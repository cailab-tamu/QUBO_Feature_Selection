import numpy as np
from scipy.sparse import csr_matrix, issparse
from joblib import Parallel, delayed
import time


def binning_method(data, method="rice"):
    """
    Number of histogram bins for `data` under the chosen rule.

    Options: "square_root", "rice", "logarithmic", "freedman-diaconis", "scott".
    Returns an int bin count.
    """
    data = np.asarray(data)
    n = data.size

    if method == "rice":
        return int(np.ceil(2 * n ** (1 / 3)))
    elif method == "logarithmic":
        return int(np.ceil(np.log2(n) + 1))
    elif method == "square_root":
        return int(np.ceil(np.sqrt(n)))
    elif method == "freedman-diaconis":
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        if iqr == 0:                       # guard: constant-ish data
            return int(np.ceil(2 * n ** (1 / 3)))
        bin_width = 2 * iqr / (n ** (1 / 3))
        return max(1, int(np.ceil((np.max(data) - np.min(data)) / bin_width)))
    elif method == "scott":
        sd = np.std(data)
        if sd == 0:
            return int(np.ceil(2 * n ** (1 / 3)))
        bin_width = (3.5 * sd) / (n ** (1 / 3))
        return max(1, int(np.ceil((np.max(data) - np.min(data)) / bin_width)))
    else:
        raise ValueError("Invalid method. Choose from: 'square_root', 'rice', "
                         "'logarithmic', 'freedman-diaconis', 'scott'.")


def _digitize_rows(matrix, bins_method, binning="equal_width"):
    """
    Extract and bin every row ONCE, returning integer bin-index rows.

    This is the key optimization: row i is densified and digitized a single time
    and then reused across all its pairings, instead of being re-extracted and
    re-binned for every partner j.

    binning:
        "equal_width" - equal-width bins over each row's min..max (matches the
                        original np.histogram2d behaviour).
        "quantile"    - equal-count (rank/quantile) bins; far more stable for
                        zero-inflated log1p expression and for MI-vs-pseudotime.
    Returns
    -------
    coded : list of (idx_array:int, nbins:int) per row
    """
    n_features = matrix.shape[0]
    coded = []
    for i in range(n_features):
        v = matrix[i, :].toarray().ravel()
        k = binning_method(v, bins_method)
        k = max(1, int(k))

        if binning == "quantile":
            # equal-count edges; handles the zero-spike gracefully
            edges = np.quantile(v, np.linspace(0, 1, k + 1))
            edges = np.unique(edges)
            if edges.size < 2:                    # constant row
                idx = np.zeros(v.size, dtype=np.int32)
                coded.append((idx, 1))
                continue
            idx = np.clip(np.digitize(v, edges[1:-1], right=False), 0,
                          edges.size - 2).astype(np.int32)
            coded.append((idx, edges.size - 1))
        else:  # equal_width
            vmin, vmax = v.min(), v.max()
            if vmax <= vmin:                       # constant row
                coded.append((np.zeros(v.size, dtype=np.int32), 1))
                continue
            edges = np.linspace(vmin, vmax, k + 1)
            idx = np.clip(np.digitize(v, edges[1:-1], right=False), 0,
                          k - 1).astype(np.int32)
            coded.append((idx, k))
    return coded


def _entropy_from_counts(counts):
    """Shannon entropy (bits) from a count array, zeros masked (no epsilon)."""
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def mutual_information_matrix(matrix, bins_method="rice", n_jobs=-1,
                              binning="equal_width", verbose=False):
    """
    Symmetric mutual-information matrix over the ROWS of `matrix`
    (features x samples), computed pairwise.

    Optimizations vs the naive version:
      * each row is densified + binned ONCE (see _digitize_rows), not per-pair
      * entropies use masked probabilities (no 1e-10 leak that inflates H)
      * bins_method is actually honoured (no hardcoded 'rice')

    Parameters
    ----------
    binning : {"equal_width", "quantile"}
        "quantile" is recommended for zero-inflated log1p expression and for
        gene-vs-pseudotime MI; "equal_width" reproduces the original behaviour.
    """
    if not issparse(matrix):
        matrix = csr_matrix(matrix)

    n_features = matrix.shape[0]
    mi_matrix = np.zeros((n_features, n_features))

    start_time = time.time()

    # --- bin every row a single time up front ---
    coded = _digitize_rows(matrix, bins_method, binning=binning)
    if verbose:
        print(f"Digitized {n_features} rows in "
              f"{time.time() - start_time:.1f}s; computing pairs...")

    def compute_pairwise_mi(i, j):
        idx_i, ki = coded[i]
        idx_j, kj = coded[j]
        # joint counts via flat 2D indexing (fast, no histogram2d re-binning)
        joint = np.bincount(idx_i * kj + idx_j, minlength=ki * kj)
        joint = joint.reshape(ki, kj)
        h_xy = _entropy_from_counts(joint)
        h_x = _entropy_from_counts(joint.sum(axis=1))
        h_y = _entropy_from_counts(joint.sum(axis=0))
        return h_x + h_y - h_xy               # MI = H(X)+H(Y)-H(X,Y)

    jobs = [(i, j) for i in range(n_features) for j in range(i + 1, n_features)]
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_pairwise_mi)(i, j) for i, j in jobs
    )

    for idx, (i, j) in enumerate(jobs):
        mi_matrix[i, j] = results[idx]
        mi_matrix[j, i] = results[idx]

    print(f"Elapsed time for MI construction: {time.time() - start_time:.4f} seconds")
    return mi_matrix



