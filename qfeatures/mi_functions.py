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
        if iqr == 0:
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


def _smallest_uint(nbins):
    """Smallest unsigned dtype that can hold bin indices 0..nbins-1."""
    if nbins <= 256:
        return np.uint8
    if nbins <= 65536:
        return np.uint16
    return np.int32


def _digitize_rows(matrix, bins_method, binning="equal_width", verbose=True):
    """
    Extract and bin every row ONCE into a compact integer code array.

    Each row is densified and digitized a single time, then reused across all of
    its pairings (instead of being re-extracted and re-binned for every partner).

    Memory: codes are stored in the smallest dtype that fits the bin count.
    For Rice on ~108k cells that is 96 bins -> uint8, so a 12k x 108k cache is
    ~1.3 GB instead of ~5.2 GB at int32.

    binning:
        "equal_width" - equal-width bins over each row's min..max.
        "quantile"    - equal-count (rank/quantile) bins; better for
                        zero-inflated log1p expression and MI-vs-pseudotime.

    Returns
    -------
    codes : ndarray (n_features, n_samples), unsigned ints
    nbins : ndarray (n_features,) int  - bins actually used per row
    """
    n_features, n_samples = matrix.shape

    # decide dtype from the nominal bin count for this sample size
    k_nominal = binning_method(np.empty(n_samples), bins_method)
    dt = _smallest_uint(max(k_nominal, 2))

    codes = np.zeros((n_features, n_samples), dtype=dt)
    nbins = np.ones(n_features, dtype=np.int64)

    t0 = time.time()
    step = max(1, n_features // 10)
    for i in range(n_features):
        v = matrix[i, :].toarray().ravel()
        k = max(1, int(binning_method(v, bins_method)))

        if binning == "quantile":
            edges = np.unique(np.quantile(v, np.linspace(0, 1, k + 1)))
            if edges.size < 2:
                nbins[i] = 1
                continue
            idx = np.clip(np.digitize(v, edges[1:-1], right=False),
                          0, edges.size - 2)
            codes[i] = idx.astype(dt, copy=False)
            nbins[i] = edges.size - 1
        else:
            vmin, vmax = v.min(), v.max()
            if vmax <= vmin:
                nbins[i] = 1
                continue
            edges = np.linspace(vmin, vmax, k + 1)
            idx = np.clip(np.digitize(v, edges[1:-1], right=False), 0, k - 1)
            codes[i] = idx.astype(dt, copy=False)
            nbins[i] = k

        if verbose and (i + 1) % step == 0:
            pct = 100 * (i + 1) / n_features
            print(f"  [binning] {pct:5.1f}%  ({i+1:,}/{n_features:,} rows, "
                  f"{time.time()-t0:.0f}s)", flush=True)

    if verbose:
        gb = codes.nbytes / 1e9
        print(f"  [binning] done in {time.time()-t0:.0f}s | cache {gb:.2f} GB "
              f"({np.dtype(dt).name}, max bins {int(nbins.max())})", flush=True)
    return codes, nbins


def _entropy_from_counts(counts):
    """Shannon entropy (bits) from a count array, zeros masked (no epsilon)."""
    total = counts.sum()
    if total == 0:
        return 0.0
    p = counts / total
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def mutual_information_matrix(matrix, bins_method="rice", n_jobs=-1,
                              binning="equal_width", backend="threading",
                              verbose=True):
    """
    Symmetric mutual-information matrix over the ROWS of `matrix`
    (features x samples).

    Notes on scale
    --------------
    Work is O(n_features^2) in pairs (12k features -> ~73M pairs), so this is
    compute-bound. Progress is reported every ~10% of rows.

    Implementation choices that matter at scale:
      * each row densified + binned ONCE, cached as uint8/uint16 (see
        _digitize_rows) rather than re-binned per pair
      * tasks are generated per-row, never materializing a ~73M-element job
        list (that alone would cost ~5 GB before any compute)
      * backend='threading' by default: numpy releases the GIL in the inner
        bincount/entropy ops, and threads SHARE the code cache instead of
        pickling it to every worker process

    Parameters
    ----------
    binning : {"equal_width", "quantile"}
        "quantile" recommended for zero-inflated log1p expression.
    backend : {"threading", "loky", None}
        joblib backend. 'threading' avoids copying the cache per worker.
    """
    if not issparse(matrix):
        matrix = csr_matrix(matrix)

    n_features = matrix.shape[0]
    n_pairs = n_features * (n_features - 1) // 2
    mi_matrix = np.zeros((n_features, n_features))

    start_time = time.time()
    if verbose:
        print(f"MI: {n_features:,} features x {matrix.shape[1]:,} samples "
              f"-> {n_pairs:,} pairs | binning={binning} ({bins_method})",
              flush=True)

    codes, nbins = _digitize_rows(matrix, bins_method, binning=binning,
                                  verbose=verbose)

    def row_mi(i):
        """MI of row i against every row j > i. Returns (i, ndarray)."""
        ci = codes[i].astype(np.int32, copy=False)
        ki = int(nbins[i])
        out = np.zeros(n_features - i - 1)
        for pos, j in enumerate(range(i + 1, n_features)):
            kj = int(nbins[j])
            joint = np.bincount(ci * kj + codes[j], minlength=ki * kj)
            joint = joint.reshape(ki, kj)
            h_xy = _entropy_from_counts(joint)
            h_x = _entropy_from_counts(joint.sum(axis=1))
            h_y = _entropy_from_counts(joint.sum(axis=0))
            out[pos] = h_x + h_y - h_xy
        return i, out

    # process rows in blocks so we can report progress ~every 10%
    rows = list(range(n_features - 1))
    n_blocks = 10
    block_size = max(1, int(np.ceil(len(rows) / n_blocks)))
    done_pairs = 0
    t_pairs = time.time()

    for b0 in range(0, len(rows), block_size):
        block = rows[b0:b0 + block_size]
        results = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(row_mi)(i) for i in block
        )
        for i, vals in results:
            mi_matrix[i, i + 1:] = vals
            mi_matrix[i + 1:, i] = vals
            done_pairs += vals.size

        if verbose:
            pct = 100 * done_pairs / n_pairs
            el = time.time() - t_pairs
            eta = el * (n_pairs - done_pairs) / max(done_pairs, 1)
            print(f"  [pairs] {pct:5.1f}%  ({done_pairs:,}/{n_pairs:,}) "
                  f"elapsed {el/60:.1f} min | ETA {eta/60:.1f} min", flush=True)

    if verbose:
        print(f"Elapsed time for MI construction: "
              f"{time.time() - start_time:.4f} seconds", flush=True)
    return mi_matrix
