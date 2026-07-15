import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, issparse


def prefilter_genes(adata, layer=None, min_cells_frac=0.01, top_n_variable=None,
                    min_dispersion=None, return_mask=False, verbose=True):
    """
    Light gene pre-filter to shrink the MI pair count BEFORE building Xy.

    The MI matrix is O(n_genes^2) in pairs, so cutting genes is the biggest
    runtime lever. This drops genes that carry little information (near-constant
    or barely-detected) with negligible loss, then optionally keeps the most
    variable ones.

    Applied to genes ONLY. The target row is never touched (append it afterward
    via mat_transform_sc as usual). Operates on the same matrix mat_transform_sc
    will read (adata.X or adata.layers[layer]) so filtering is consistent.

    Parameters
    ----------
    adata : AnnData  (cells x genes)
    layer : str or None
        Which matrix to compute statistics on (default adata.X).
    min_cells_frac : float
        Drop genes detected (nonzero) in fewer than this FRACTION of cells.
        0.01 = expressed in <1% of cells. Set 0 to disable.
    top_n_variable : int or None
        If set, additionally keep only the top-N genes by variance (computed on
        the log1p values). E.g. 3000. If None, keep all that pass the detection
        and dispersion filters.
    min_dispersion : float or None
        If set, drop genes whose variance/mean (index of dispersion) is below
        this. Cheap way to remove flat, uninformative genes. If None, skip.
    return_mask : bool
        If True, also return the boolean gene mask (length n_genes of the input).

    Returns
    -------
    adata_sub : AnnData
        View/copy restricted to the kept genes. Pass this straight into
        mat_transform_sc.
    mask : np.ndarray (only if return_mask=True)
        Boolean over the ORIGINAL genes indicating which were kept.
    """
    X = adata.layers[layer] if layer is not None else adata.X   # cells x genes
    if not issparse(X):
        X = csr_matrix(X)
    X = X.tocsc()                       # column ops (per-gene) are fast in CSC
    n_cells, n_genes = X.shape

    keep = np.ones(n_genes, dtype=bool)

    # 1) detection filter: nonzero in at least min_cells_frac of cells
    if min_cells_frac and min_cells_frac > 0:
        n_nonzero = np.diff(X.indptr)               # nonzeros per gene (column)
        det = n_nonzero >= (min_cells_frac * n_cells)
        keep &= det

    # per-gene mean / variance via E[x^2] - E[x]^2 (sparse-safe, no densify)
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = np.maximum(mean_sq - mean ** 2, 0.0)

    # 2) dispersion filter (variance / mean), guards divide-by-zero
    if min_dispersion is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            disp = np.where(mean > 0, var / mean, 0.0)
        keep &= disp >= min_dispersion

    # 3) optional top-N by variance among the survivors
    if top_n_variable is not None and keep.sum() > top_n_variable:
        surv = np.flatnonzero(keep)
        order = surv[np.argsort(var[surv])[::-1]]   # high variance first
        top = order[:top_n_variable]
        new_keep = np.zeros(n_genes, dtype=bool)
        new_keep[top] = True
        keep = new_keep

    if verbose:
        print(f"[prefilter_genes] {n_genes} -> {int(keep.sum())} genes kept "
              f"({n_genes - int(keep.sum())} dropped). "
              f"MI pairs: {n_genes*(n_genes-1)//2:,} -> "
              f"{int(keep.sum())*(int(keep.sum())-1)//2:,}")

    adata_sub = adata[:, keep].copy()
    if return_mask:
        return adata_sub, keep
    return adata_sub


def compute_pearson_residual(matrix, theta=100):
    """
    Compute clipped analytic Pearson residuals for a sparse matrix with rows as
    features and columns as cells.
    Reference: https://doi.org/10.1101/2020.12.01.405886

    WARNING - MEMORY: Pearson residuals are DENSE by construction. Every zero
    count maps to a nonzero residual, so the output has no zeros to drop. The
    intermediate `expected` outer product alone is a full features x cells dense
    array (~10 GB in float64 for 12k x 108k), and several dense copies are made
    during the computation (peak tens of GB). Do NOT call this on a full
    12k-gene x 108k-cell matrix. Only use it on a gene SUBSET (e.g. 1-2k HVGs):
        hvg = adata.var["highly_variable"].to_numpy()
        res = compute_pearson_residual(adata[:, hvg].X.T)   # features x cells
    For the full matrix, use mat_transform_sc(..., method="log1p") instead.

    Args:
        matrix (scipy.sparse.csr_matrix or numpy.ndarray): rows = features,
            columns = cells.
        theta (float): NB overdispersion parameter.

    Returns:
        scipy.sparse.csr_matrix: clipped Pearson residuals (dense-valued).
    """
    if not issparse(matrix):
        matrix = csr_matrix(matrix)

    row_sums = matrix.sum(axis=1)
    col_sums = matrix.sum(axis=0)
    total_sum = matrix.sum()

    row_sums_dense = np.array(row_sums).flatten()
    col_sums_dense = np.array(col_sums).flatten()
    expected = (row_sums_dense[:, None] * col_sums_dense[None, :]) / total_sum

    expected_sparse = csr_matrix(expected)
    std_dev = np.sqrt(expected_sparse + (expected_sparse.power(2)) / theta)

    residuals = (matrix - expected_sparse).multiply(std_dev.power(-1))
    residuals.data[np.isnan(residuals.data)] = 0
    residuals.data[np.isinf(residuals.data)] = 0

    num_cells = matrix.shape[1]
    clip_value = np.sqrt(num_cells)
    residuals.data = np.clip(residuals.data, -clip_value, clip_value)
    return residuals


def _calibrate_mi_rate(n_cells, n_bins=96, n_trials=200):
    """
    Measure how many MI pairs/sec this machine does, by timing the same inner
    loop mi_functions uses (bincount on binned rows + entropies). Single
    threaded; the caller scales by n_jobs.
    """
    import time
    rng = np.random.default_rng(0)
    a = rng.integers(0, n_bins, n_cells).astype(np.int32)
    b = rng.integers(0, n_bins, n_cells).astype(np.uint8)

    def _ent(c):
        t = c.sum()
        if t == 0:
            return 0.0
        p = c / t
        p = p[p > 0]
        return float(-np.sum(p * np.log2(p)))

    t0 = time.time()
    for _ in range(n_trials):
        joint = np.bincount(a * n_bins + b, minlength=n_bins * n_bins)
        joint = joint.reshape(n_bins, n_bins)
        _ent(joint); _ent(joint.sum(axis=1)); _ent(joint.sum(axis=0))
    el = time.time() - t0
    return n_trials / el if el > 0 else float("inf")


def sweep_prefilter(adata, layer=None, n_jobs=1, target_minutes=30,
                    calibrate=True, detect_fracs=(0.001, 0.005, 0.01, 0.05, 0.1),
                    disp_thresholds=(0.1, 0.2, 0.3, 0.5, 1.0, 2.0),
                    top_ns=(10000, 7000, 5000, 3000, 2000)):
    """
    Diagnose gene-filter options before committing to an MI run.

    Prints, for the matrix mat_transform_sc will read:
      * detection and dispersion distributions
      * genes kept / MI pairs / ESTIMATED MINUTES for each candidate threshold
      * a recommendation that fits `target_minutes`

    The runtime estimate is calibrated by timing this machine on the same inner
    loop mutual_information_matrix uses, so it reflects your hardware rather
    than a guess. It is still an estimate: real runs vary with thread
    contention and the per-row binning pass.

    Parameters
    ----------
    adata : AnnData (cells x genes)
    layer : str or None
        Matrix to analyse (default adata.X). Should be the sparse log1p data.
    n_jobs : int
        Threads you plan to pass to mutual_information_matrix. Used only to
        scale the time estimate. Pass the real number (-1 -> os.cpu_count()).
    target_minutes : float
        Runtime budget used to pick the recommendation.
    calibrate : bool
        If False, skip the micro-benchmark and report pairs only.

    Returns
    -------
    dict with keys 'mean', 'var', 'dispersion', 'detect_frac', 'pairs_per_sec'
    so you can do your own analysis.
    """
    import os

    X = adata.layers[layer] if layer is not None else adata.X
    if not issparse(X):
        X = csr_matrix(X)
    X = X.tocsc()
    n_cells, n_genes = X.shape

    # per-gene stats, sparse-safe
    n_nonzero = np.diff(X.indptr)
    detect_frac = n_nonzero / n_cells
    mean = np.asarray(X.mean(axis=0)).ravel()
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    var = np.maximum(mean_sq - mean ** 2, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        disp = np.where(mean > 0, var / mean, 0.0)

    if n_jobs is None or n_jobs < 0:
        n_jobs_eff = os.cpu_count() or 1
    else:
        n_jobs_eff = max(1, n_jobs)

    pps = None
    if calibrate:
        rate1 = _calibrate_mi_rate(n_cells)
        pps = rate1 * n_jobs_eff          # optimistic: perfect thread scaling
        print(f"[sweep] calibrated ~{rate1:,.0f} pairs/s/thread x {n_jobs_eff} "
              f"threads -> ~{pps:,.0f} pairs/s (optimistic)")

    def fmt_time(n_pairs):
        if pps is None:
            return ""
        mins = n_pairs / pps / 60
        if mins < 90:
            return f"{mins:>8.1f} min"
        return f"{mins/60:>8.1f} hr "

    print(f"\n=== {n_genes:,} genes x {n_cells:,} cells ===")
    print(f"all genes -> {n_genes*(n_genes-1)//2:,} pairs {fmt_time(n_genes*(n_genes-1)//2)}")

    print("\ndetection (fraction of cells with nonzero):")
    for p in (10, 25, 50, 75, 90):
        print(f"  p{p:<2}: {np.percentile(detect_frac, p):.3f}")

    print("\ndispersion (var/mean):")
    dp = {p: np.percentile(disp, p) for p in (10, 25, 50, 75, 90)}
    for p, v in dp.items():
        print(f"  p{p:<2}: {v:.3f}")

    # is dispersion actually discriminative here?
    spread = dp[90] / dp[10] if dp[10] > 0 else np.inf
    disp_useless = spread < 2.5
    if disp_useless:
        print(f"  !! p90/p10 = {spread:.2f} -> dispersion is NOT discriminative.")
        print("     On log1p data var/mean collapses toward a common value, so any")
        print("     threshold here is arbitrary (tiny changes swing the gene count")
        print("     wildly). Prefer detection + top_n_variable.")

    print("\nmin_cells_frac (detection filter):")
    for f in detect_fracs:
        n = int((detect_frac >= f).sum())
        pr = n * (n - 1) // 2
        note = "  <- drops rare-population markers" if f >= 0.05 else ""
        print(f"  {f:<6} -> {n:>6,} genes, {pr:>12,} pairs {fmt_time(pr)}{note}")

    print("\nmin_dispersion:")
    for thr in disp_thresholds:
        n = int((disp >= thr).sum())
        pr = n * (n - 1) // 2
        print(f"  {thr:<6} -> {n:>6,} genes, {pr:>12,} pairs {fmt_time(pr)}")

    print("\ntop_n_variable (rank by variance):")
    for n in top_ns:
        if n > n_genes:
            continue
        pr = n * (n - 1) // 2
        print(f"  {n:<6} -> {n:>6,} genes, {pr:>12,} pairs {fmt_time(pr)}")

    # ---- recommendation ----
    print("\n=== recommendation ===")
    if pps is None:
        print("  (no calibration; re-run with calibrate=True for time estimates)")
    else:
        budget_pairs = pps * target_minutes * 60
        n_ok = int((1 + np.sqrt(1 + 8 * budget_pairs)) / 2)   # solve n(n-1)/2<=B
        n_ok = min(n_ok, n_genes)
        print(f"  budget {target_minutes:g} min -> about {n_ok:,} genes "
              f"({n_ok*(n_ok-1)//2:,} pairs)")
        cand = [n for n in top_ns if n <= n_ok]
        pick = max(cand) if cand else min(top_ns)
        base = "min_cells_frac=0.01"
        if n_ok >= n_genes:
            print(f"  all {n_genes:,} genes fit the budget; "
                  f"prefilter_genes(adata, {base}) is enough.")
        else:
            print(f"  prefilter_genes(adata, {base}, top_n_variable={pick})")
        if disp_useless:
            print("  (omit min_dispersion: not discriminative on this data)")
        print("\n  NOTE: top_n_variable ranks by variance, which is a crude HVG "
              "selection.\n  It is a choice about your feature space, not a "
              "principled cutoff, so\n  record whatever you use in your methods. "
              "scanpy's highly_variable_genes\n  does the same job better "
              "(it bins by mean first).")

    return {"mean": mean, "var": var, "dispersion": disp,
            "detect_frac": detect_frac, "pairs_per_sec": pps}


def _mem_report(name, Xsp):
    """Print shape / nnz / approximate CSR memory footprint."""
    total = Xsp.data.nbytes + Xsp.indices.nbytes + Xsp.indptr.nbytes
    print(f"{name}: shape={Xsp.shape}  nnz={Xsp.nnz:,}  "
          f"~{total / 1e9:.2f} GB  dtype={Xsp.data.dtype}")


def verify_expression(adata, layer=None):
    """
    Sanity-check that the matrix mat_transform_sc will read is sparse
    log1p-normalized data (and NOT scaled/centered data, which densifies and
    breaks the sparse assumption).

    Prints a short report and returns True if it looks like log1p, False if it
    looks scaled (has negatives) or is dense.
    """
    X = adata.layers[layer] if layer is not None else adata.X
    sparse = issparse(X)
    xmin = float(X.min())
    xmax = float(X.max())
    where = f"adata.layers['{layer}']" if layer is not None else "adata.X"
    looks_log1p = sparse and xmin >= 0.0
    print(f"{where}: sparse={sparse}  min={xmin:.3f}  max={xmax:.3f}")
    if not sparse:
        print("  -> DENSE. Likely scaled; keep a sparse lognorm layer and pass layer=.")
    elif xmin < 0.0:
        print("  -> has negatives -> looks SCALED/centered, not log1p. Use a lognorm layer.")
    else:
        print("  -> looks like log1p + libsize (non-negative, sparse). Good to use.")
    return looks_log1p


def mat_transform_sc(adata, obs_key, layer=None, method="log1p", theta=100,
                     dtype=np.float32, drop_nan_target=True):
    """
    Build the (n_features + 1) x n_cells feature-by-cell matrix for QUBO /
    mutual-information feature selection: gene rows stacked on top of a single
    target row taken from adata.obs[obs_key] (e.g. a pseudotime column).

    Stays sparse end to end and never densifies the full matrix.

    Parameters
    ----------
    adata : AnnData
    obs_key : str
        Column in adata.obs to use as the target (last row of Xy).
    layer : str or None
        Use adata.layers[layer] instead of adata.X. Point this at your
        log1p-normalized layer if adata.X has been scaled/centered (scaling
        densifies, which defeats the purpose).
    method : {"log1p", "pearson"}
        "log1p" (default): use adata.X (or the given layer) AS IS. This does NOT
            re-normalize - it assumes X already holds log1p + libsize-normalized
            values. Run verify_expression(adata) first if unsure. Cheap and
            sparse (~1 GB for 12k x 108k at typical scRNA-seq density).
        "pearson": analytic Pearson residuals. DENSE and expensive; only sane on
            a gene subset. See compute_pearson_residual's warning.
    dtype : np.dtype
        float32 halves memory vs float64 and is plenty for binned MI.
    drop_nan_target : bool
        Drop cells whose target is NaN (e.g. unreachable pseudotime) so the MI
        step never sees NaNs. Columns are dropped from Xy.

    Returns
    -------
    Xy : scipy.sparse.csr_matrix, shape (n_genes + 1, n_cells_kept)
        Rows 0..n_genes-1 are genes; the last row is the target.
    """
    X = adata.layers[layer] if layer is not None else adata.X  # cells x genes

    y = np.asarray(adata.obs[obs_key].values, dtype=np.float64)

    # drop cells with NaN target (keeps MI happy; columns of the final Xy)
    if drop_nan_target and np.isnan(y).any():
        keep = np.isfinite(y)
        print(f"[mat_transform_sc] dropping {int((~keep).sum())} cell(s) with "
              f"NaN '{obs_key}' ({int(keep.sum())} kept)")
        X = X[keep]
        y = y[keep]

    if method == "log1p":
        if not issparse(X):
            X = csr_matrix(X)
        # genes x cells, sparse, never densified
        Xf = X.T.tocsr().astype(dtype)
    elif method == "pearson":
        # X.T is features x cells; DENSE inside - only do this on a gene subset
        Xf = csr_matrix(compute_pearson_residual(X.T, theta=theta)).astype(dtype)
    else:
        raise ValueError("method must be 'log1p' or 'pearson'")

    # append the target as one sparse row (only n_cells nonzeros added)
    y_row = csr_matrix(y.reshape(1, -1).astype(dtype))
    Xy = sp.vstack([Xf, y_row], format="csr")

    _mem_report("features", Xf)
    _mem_report("Xy (final)", Xy)
    return Xy


def mat_transform(X, y, method="log1p", theta=100, dtype=np.float32):
    """
    Same as mat_transform_sc but for a raw expression matrix X (cells x genes)
    and a target vector y, without an AnnData wrapper. Stays sparse.

    Returns
    -------
    Xy : scipy.sparse.csr_matrix, shape (n_genes + 1, n_cells)
    """
    if method == "log1p":
        if not issparse(X):
            X = csr_matrix(X)
        Xf = X.T.tocsr().astype(dtype)          # genes x cells
    elif method == "pearson":
        Xf = csr_matrix(compute_pearson_residual(X.T, theta=theta)).astype(dtype)
    else:
        raise ValueError("method must be 'log1p' or 'pearson'")

    y_row = csr_matrix(np.asarray(y, dtype=dtype).reshape(1, -1))
    Xy = sp.vstack([Xf, y_row], format="csr")

    _mem_report("features", Xf)
    _mem_report("Xy (final)", Xy)
    return Xy
