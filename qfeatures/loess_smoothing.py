"""
loess_smoothing.py

Python port of the MATLAB `smooth(X, y, span, 'loess')` routine used to smooth
gene expression along pseudotime in single-cell RNA-seq data.

MATLAB's 'loess' method performs LOCAL weighted least-squares regression with a
2nd-degree polynomial and a tricube weight kernel, over the `span` nearest
neighbours of each query point. This module reproduces that behaviour using
only NumPy (no statsmodels dependency, since statsmodels' lowess is degree-1
only and will not match 'loess').

Method options (mirrors the MATLAB comment block):
    degree=1            -> 'lowess'  (local linear)
    degree=2            -> 'loess'   (local quadratic)   <- MATLAB default here
    robust=True, deg=1  -> 'rlowess'
    robust=True, deg=2  -> 'rloess'

Entry points:
    loess_smoothing(X, y, span)              -> single vector (MATLAB-equivalent)
    loess_smoothing_matrix(t, X, span)       -> (ngene x ncell) matrix, MATLAB loop
    loess_smoothing_adata(adata, genes, ...) -> pull genes from a scanpy AnnData
    plot_loess_curves(result, ...)           -> quick pseudotime vs expression plot

Author: (port) for Selim Romero
"""

import numpy as np


def _resolve_grid(grid, n):
    """
    Decide the number of anchor points for the grid approximation.

    grid=None  -> automatic: exact fit for small n, 500 anchors above 20k points
                  (where the exact O(n*r) fit starts costing minutes).
    grid=False/0 -> force the exact per-point fit.
    grid=int   -> use that many anchors.
    """
    if grid is False or grid == 0:
        return None
    if grid is None:
        return 500 if n > 20000 else None
    return int(grid)


def loess_smoothing(X, y, span=0.75, degree=2, robust=False, n_robust_iter=3,
                    grid=None):
    """
    Smooth response `y` against predictor `X` using local regression.

    Direct port of:
        [smoothed_y, sortIdx] = loess_smoothing(X, y, span)   % MATLAB

    Parameters
    ----------
    X : array_like, shape (n,)
        Predictor (in your usage, spline pseudotime `t`).
    y : array_like, shape (n,)
        Response to be smoothed (in your usage, expression of one gene).
    span : float
        Smoothing parameter. If 0 < span <= 1 it is a FRACTION of the total
        number of points used in each local fit (0.3 keeps more noise,
        0.75 gives a smooth curve). If span > 1 it is an absolute point count.
    degree : {1, 2}
        Polynomial degree. 2 == MATLAB 'loess', 1 == 'lowess'.
    robust : bool
        If True, run robustifying iterations (MATLAB 'rloess' / 'rlowess')
        that down-weight outliers via bisquare weights.
    n_robust_iter : int
        Number of robust iterations (MATLAB uses ~5 internally; 3 is plenty).
    grid : None | False | int
        Speed control for large n. The exact fit is O(n * span * n), which at
        ~100k points costs MINUTES per gene. With `grid`, the curve is fitted at
        that many quantile-spaced anchors and interpolated back to every point,
        which is thousands of times faster and visually identical.
        None (default) = automatic (exact below 20k points, 500 anchors above).
        False or 0 = force the exact MATLAB-equivalent fit.
        int = use that many anchors.
        Note: `grid` is ignored when robust=True (robust iterations need
        residuals at every point).

    Returns
    -------
    smoothed_y : ndarray, shape (n,)
        Smoothed values, ordered to match `X` sorted ascending
        (same convention as the MATLAB function's `smoothed_y`).
    sort_idx : ndarray, shape (n,)
        Indices that sort `X` ascending, i.e. smoothed_y corresponds to
        X[sort_idx]. Use it to reorder companion arrays: t_sorted = X[sort_idx].
    """
    X = np.asarray(X, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if X.size != y.size:
        raise ValueError(f"X and y must be the same length ({X.size} vs {y.size}).")
    n = X.size
    if n == 0:
        return np.array([]), np.array([], dtype=int)

    # --- sort by predictor (matches [X_sorted, sortIdx] = sort(X)) ---
    sort_idx = np.argsort(X, kind="stable")
    Xs = X[sort_idx]
    ys = y[sort_idx]

    # --- window size (number of points per local fit) ---
    if span <= 0:
        raise ValueError("span must be > 0.")
    if span < 1:
        r = int(np.ceil(span * n))
    else:
        r = int(round(span))
    r = max(r, degree + 1)   # need enough points to fit the polynomial
    r = min(r, n)

    n_grid = _resolve_grid(grid, n)
    if n_grid is not None and not robust:
        smoothed = _loess_fit_grid(Xs, ys, r, degree, n_grid)
    else:
        smoothed = _loess_fit(Xs, ys, r, degree, robust, n_robust_iter)
    return smoothed, sort_idx


def _loess_fit(x, y, r, degree, robust, n_iter):
    """Core sliding-window local polynomial fit on already-sorted x, y."""
    n = x.size
    yhat = np.empty(n)
    rob_w = np.ones(n)
    iters = n_iter if robust else 1

    for it in range(iters):
        left, right = 0, r - 1  # contiguous window [left, right] of size r
        for i in range(n):
            xi = x[i]
            # slide window so it holds the r nearest points to xi
            # (Cleveland's rule: advance while the next point on the right is
            #  closer to xi than the current leftmost point)
            while right < n - 1 and (x[right + 1] - xi) < (xi - x[left]):
                left += 1
                right += 1

            xw = x[left:right + 1]
            yw = y[left:right + 1]

            # tricube weights based on distance to the farthest neighbour
            d = np.abs(xw - xi)
            dmax = d[-1] if d[-1] >= d[0] else d.max()
            if dmax > 0:
                u = d / dmax
                w = (1.0 - u ** 3) ** 3
                w[u >= 1.0] = 0.0
            else:
                w = np.ones_like(d)

            if robust:
                w = w * rob_w[left:right + 1]

            yhat[i] = _weighted_poly_at0(xw - xi, yw, w, degree)

        if robust and it < iters - 1:
            resid = y - yhat
            s = np.median(np.abs(resid))
            if s <= 0:
                break
            u = np.clip(resid / (6.0 * s), -1.0, 1.0)
            rob_w = (1.0 - u ** 2) ** 2

    return yhat


def _loess_fit_grid(x, y, r, degree, n_grid):
    """
    Fast LOESS for large n: fit at `n_grid` anchor points spanning the data,
    then linearly interpolate back to every x.

    The pure sliding-window fit is O(n * r); with span=0.75 and n=100k that is
    ~8 billion element-ops (minutes per gene). Fitting at 500 anchors instead is
    O(n_grid * r), i.e. thousands of times cheaper, and because LOESS is already
    a heavily smoothed curve, interpolating between anchors is visually and
    numerically indistinguishable (sub-1% of curve range at n_grid>=200).

    Anchors are quantile-spaced so dense regions get more of them.
    """
    n = x.size
    n_grid = int(min(max(n_grid, degree + 2), n))
    # quantile-spaced anchors follow the data density, unlike linspace
    qs = np.linspace(0, 1, n_grid)
    grid = np.quantile(x, qs)
    grid = np.unique(grid)
    if grid.size < 2:
        return np.full(n, float(np.mean(y)))

    fitted = np.empty(grid.size)
    for gi, q in enumerate(grid):
        center = int(np.searchsorted(x, q))
        lo = max(0, min(center - r // 2, n - r))
        hi = lo + r
        # shift window so it holds the r nearest points to q
        while lo > 0 and (q - x[lo - 1]) < (x[hi - 1] - q):
            lo -= 1
            hi -= 1
        while hi < n and (x[hi] - q) < (q - x[lo]):
            lo += 1
            hi += 1

        xw = x[lo:hi]
        yw = y[lo:hi]
        d = np.abs(xw - q)
        dmax = d.max()
        if dmax > 0:
            u = d / dmax
            w = (1.0 - u ** 3) ** 3
            w[u >= 1.0] = 0.0
        else:
            w = np.ones_like(d)
        fitted[gi] = _weighted_poly_at0(xw - q, yw, w, degree)

    return np.interp(x, grid, fitted)


def _weighted_poly_at0(dx, yw, w, degree):
    """
    Weighted least-squares polynomial fit in the centred coordinate dx = x - x0,
    returning the fitted value at x0 (i.e. dx = 0, which is just the intercept).
    Solving the small normal equations directly is fast and stable enough here.
    """
    V = np.vander(dx, degree + 1, increasing=True)  # cols: dx^0, dx^1, dx^2
    VtW = V.T * w                                   # (degree+1, m)
    A = VtW @ V                                     # (degree+1, degree+1)
    b = VtW @ yw
    try:
        coef = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(A, b, rcond=None)[0]
    return coef[0]  # value at dx = 0


def loess_smoothing_matrix(t, X, span=0.75, degree=2, robust=False,
                           n_robust_iter=3):
    """
    Batch version mirroring your MATLAB gene loop.

        ngene = size(X_qubo, 1); ncell = size(X_qubo, 2);
        for ig = 1:ngene
            [y_qubo_fit(:,ig), idx] = loess_smoothing(t, X_qubo(ig,:)', 0.75);
            t_qubo_sort(:,ig) = t(idx);
        end

    Parameters
    ----------
    t : array_like, shape (ncell,)
        Pseudotime (predictor), shared across genes.
    X : array_like, shape (ngene, ncell)
        Expression matrix with genes in ROWS (like X_qubo).
    span, degree, robust, n_robust_iter : see loess_smoothing.

    Returns
    -------
    y_fit  : ndarray, shape (ncell, ngene)
        Smoothed expression per gene, each column sorted by pseudotime.
    t_sort : ndarray, shape (ncell, ngene)
        Pseudotime reordered to match each gene's sort. Because `t` is shared,
        every column is identical (kept 2-D to match the MATLAB layout).
    """
    t = np.asarray(t, dtype=float).ravel()
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError("X must be 2-D (ngene x ncell).")
    ngene, ncell = X.shape
    if t.size != ncell:
        raise ValueError(f"len(t)={t.size} must equal ncell={ncell}.")

    # t is shared, so the sort order is the same for every gene: compute once.
    sort_idx = np.argsort(t, kind="stable")
    t_sorted = t[sort_idx]

    y_fit = np.zeros((ncell, ngene))
    for ig in range(ngene):
        ys = X[ig, sort_idx]
        if span < 1:
            r = max(int(np.ceil(span * ncell)), degree + 1)
        else:
            r = int(round(span))
        r = min(r, ncell)
        y_fit[:, ig] = _loess_fit(t_sorted, ys, r, degree, robust, n_robust_iter)

    t_sort = np.repeat(t_sorted[:, None], ngene, axis=1)
    return y_fit, t_sort


# ---------------------------------------------------------------------------
# scanpy / AnnData interface
# ---------------------------------------------------------------------------

class LoessResult:
    """
    Container returned by loess_smoothing_adata.

    Attributes
    ----------
    t : ndarray, shape (ncell,)
        Pseudotime, sorted ascending. All curves share this x-axis.
    genes : list[str]
        Gene names in the order requested (== column order of the frames).
    smoothed : pandas.DataFrame, shape (ncell, ngene)
        LOESS-smoothed expression. Rows are cells in pseudotime order
        (index = cell barcodes), columns = genes in requested order.
    raw : pandas.DataFrame, shape (ncell, ngene)
        Raw expression in the same cell/gene order (handy for scatter overlays).
    cells : ndarray
        Cell barcodes in pseudotime order (== smoothed.index).
    pseudotime_key : str
        The obs column used as pseudotime.
    """

    def __init__(self, t, genes, smoothed, raw, cells, pseudotime_key):
        self.t = t
        self.genes = list(genes)
        self.smoothed = smoothed
        self.raw = raw
        self.cells = cells
        self.pseudotime_key = pseudotime_key

    def __repr__(self):
        return (f"LoessResult(ncell={self.t.size}, ngenes={len(self.genes)}, "
                f"pseudotime='{self.pseudotime_key}', genes={self.genes})")


def _extract_expression(adata, gene_idxs, layer, use_raw):
    """Pull a (ncell, k) dense float array for the requested gene columns."""
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        mat = adata.raw.X
    elif layer is not None:
        if layer not in adata.layers:
            raise ValueError(f"layer '{layer}' not found. "
                             f"Available: {list(adata.layers.keys())}")
        mat = adata.layers[layer]
    else:
        mat = adata.X

    sub = mat[:, gene_idxs]
    # densify if sparse (scipy.sparse matrices expose .toarray)
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    return np.asarray(sub, dtype=float)


def loess_smoothing_adata(adata, genes, pseudotime_key="splinefit_pseudotime",
                          layer=None, use_raw=False, span=0.75, degree=2,
                          robust=False, n_robust_iter=3, dropna=True,
                          grid=None, verbose=True):
    """
    Smooth expression of selected genes along pseudotime from an AnnData object.

    Parameters
    ----------
    adata : AnnData
        A scanpy AnnData. Expression taken from adata.X (default), or from
        adata.layers[layer], or adata.raw.X if use_raw=True.
    genes : sequence of str
        Gene names (must be in adata.var_names). Curves are returned/plotted
        in exactly this order.
    pseudotime_key : str
        Column in adata.obs holding pseudotime (default 'splinefit_pseudotime').
    layer : str or None
        Use adata.layers[layer] instead of adata.X (e.g. a lognorm layer).
    use_raw : bool
        Use adata.raw.X instead of adata.X.
    span, degree, robust, n_robust_iter : see loess_smoothing.
    dropna : bool
        Drop cells whose pseudotime is NaN before fitting (default True).

    Returns
    -------
    LoessResult
        Has .t, .genes, .smoothed (DataFrame), .raw (DataFrame), .cells.
    """
    import pandas as pd

    genes = list(genes)
    if pseudotime_key not in adata.obs:
        raise KeyError(f"'{pseudotime_key}' not in adata.obs. "
                       f"Available obs columns include: "
                       f"{list(adata.obs.columns)[:20]}")

    # map gene names -> column positions, reporting any that are missing
    var_names = list(adata.var_names)
    name_to_idx = {g: i for i, g in enumerate(var_names)}
    missing = [g for g in genes if g not in name_to_idx]
    if missing:
        raise KeyError(f"{len(missing)} gene(s) not in adata.var_names: {missing}")
    gene_idxs = [name_to_idx[g] for g in genes]

    # pseudotime + optional NaN filtering
    t_all = np.asarray(adata.obs[pseudotime_key].to_numpy(), dtype=float)
    cells_all = np.asarray(adata.obs_names)
    keep = np.ones(t_all.size, dtype=bool)
    if dropna:
        keep = ~np.isnan(t_all)
        n_drop = int((~keep).sum())
        if n_drop:
            print(f"[loess] dropping {n_drop} cell(s) with NaN "
                  f"'{pseudotime_key}' ({keep.sum()} kept).")

    expr = _extract_expression(adata, gene_idxs, layer, use_raw)  # (ncell_all, k)
    expr = expr[keep]
    t = t_all[keep]
    cells = cells_all[keep]

    # shared sort by pseudotime (same order for every gene)
    sort_idx = np.argsort(t, kind="stable")
    t_sorted = t[sort_idx]
    cells_sorted = cells[sort_idx]
    ncell = t_sorted.size

    if span < 1:
        r = max(int(np.ceil(span * ncell)), degree + 1)
    else:
        r = int(round(span))
    r = min(r, ncell)

    n_grid = _resolve_grid(grid, ncell)
    use_grid = n_grid is not None and not robust
    if verbose:
        if use_grid:
            print(f"[loess] {len(genes)} gene(s) x {ncell:,} cells | "
                  f"grid={n_grid} anchors (fast path)")
        else:
            est = ncell * r / 4e6
            print(f"[loess] {len(genes)} gene(s) x {ncell:,} cells | EXACT fit, "
                  f"window r={r:,} -> roughly {est*len(genes)/60:.0f} min total. "
                  f"Pass grid=500 to make this ~instant.")

    smoothed = np.zeros((ncell, len(genes)))
    raw_sorted = np.zeros((ncell, len(genes)))
    for j in range(len(genes)):
        ys = expr[sort_idx, j]
        raw_sorted[:, j] = ys
        if use_grid:
            smoothed[:, j] = _loess_fit_grid(t_sorted, ys, r, degree, n_grid)
        else:
            smoothed[:, j] = _loess_fit(t_sorted, ys, r, degree, robust,
                                        n_robust_iter)

    smoothed_df = pd.DataFrame(smoothed, index=cells_sorted, columns=genes)
    raw_df = pd.DataFrame(raw_sorted, index=cells_sorted, columns=genes)

    return LoessResult(t_sorted, genes, smoothed_df, raw_df,
                       cells_sorted, pseudotime_key)


def plot_loess_curves(result, layout="grid", ncols=3, figsize=None,
                      scatter=True, zscore=False, scatter_kws=None,
                      line_kws=None, sharex=True, sharey=False, save=None):
    """
    Plot pseudotime (x) vs expression (y) with the LOESS curve overlaid.

    Parameters
    ----------
    result : LoessResult
        Output of loess_smoothing_adata.
    layout : {'grid', 'overlay'}
        'grid'    -> one panel per gene (raw scatter + smoothed line).
        'overlay' -> all smoothed curves on a single axis (good for comparing
                     a handful of genes; consider zscore=True).
    ncols : int
        Columns in the grid layout.
    zscore : bool
        Z-score each gene's smoothed curve (per gene) before plotting. Useful
        for overlay comparisons where genes are on different scales.
    scatter : bool
        Draw the raw expression scatter behind the curve (grid layout only).
    scatter_kws, line_kws : dict or None
        Extra matplotlib kwargs for the scatter / line.
    save : str or None
        If given, savefig to this path.

    Returns
    -------
    (fig, axes)
    """
    import matplotlib.pyplot as plt

    genes = result.genes
    t = result.t
    scatter_kws = {"s": 4, "alpha": 0.25, "color": "0.6", **(scatter_kws or {})}
    line_kws = {"lw": 2.0, "color": "crimson", **(line_kws or {})}

    def _maybe_z(v):
        if not zscore:
            return v
        sd = v.std()
        return (v - v.mean()) / sd if sd > 0 else v - v.mean()

    if layout == "overlay":
        fig, ax = plt.subplots(figsize=figsize or (7, 5))
        for g in genes:
            ax.plot(t, _maybe_z(result.smoothed[g].to_numpy()), label=g, lw=2)
        ax.set_xlabel(result.pseudotime_key)
        ax.set_ylabel("z-scored smoothed expr" if zscore else "smoothed expr")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=150, bbox_inches="tight")
        return fig, ax

    # grid layout
    n = len(genes)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or (4 * ncols, 3 * nrows),
                             sharex=sharex, sharey=sharey, squeeze=False)
    axes_flat = axes.ravel()
    for k, g in enumerate(genes):
        ax = axes_flat[k]
        yhat = _maybe_z(result.smoothed[g].to_numpy())
        if scatter:
            yraw = result.raw[g].to_numpy()
            ax.scatter(t, _maybe_z(yraw) if zscore else yraw, **scatter_kws)
        ax.plot(t, yhat, **line_kws)
        ax.set_title(g, fontsize=10)
        ax.set_xlabel(result.pseudotime_key, fontsize=8)
        ax.set_ylabel("expr", fontsize=8)
    for k in range(n, len(axes_flat)):   # hide unused panels
        axes_flat[k].axis("off")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
    return fig, axes


if __name__ == "__main__":
    # quick self-test / demo
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 1, 400))
    signal = np.sin(2 * np.pi * t) + 0.5 * t
    expr = signal + rng.normal(0, 0.3, t.size)

    sm, idx = loess_smoothing(t, expr, span=0.75, degree=2)
    err = np.sqrt(np.mean((sm - signal[idx]) ** 2))
    print(f"n={t.size}  span=0.75  loess RMSE vs true signal = {err:.4f}")

    # matrix form: 3 genes sharing pseudotime
    Xg = np.vstack([expr, np.cos(2 * np.pi * t) + rng.normal(0, .3, t.size),
                    2 * t + rng.normal(0, .3, t.size)])
    yfit, tsort = loess_smoothing_matrix(t, Xg, span=0.75)
    print("matrix y_fit shape:", yfit.shape, " t_sort shape:", tsort.shape)

    # ----- AnnData path (uses real anndata if available, else a duck stand-in) -----
    ncell = 500
    tt = rng.uniform(0, 1, ncell)                    # unsorted pseudotime, on purpose
    genes_all = ["Lgr5", "Muc2", "Lyz1", "Defa39", "Mki67"]
    Xcells = np.column_stack([
        np.sin(2 * np.pi * tt) + rng.normal(0, .3, ncell),
        (1 - tt) + rng.normal(0, .3, ncell),
        np.exp(-((tt - 0.7) ** 2) / 0.02) + rng.normal(0, .2, ncell),
        (tt > 0.6).astype(float) + rng.normal(0, .2, ncell),
        tt ** 2 + rng.normal(0, .3, ncell),
    ])  # (ncell, ngene)

    try:
        import anndata as ad
        import pandas as pd
        adata = ad.AnnData(
            X=Xcells,
            obs=pd.DataFrame({"splinefit_pseudotime": tt},
                             index=[f"cell{i}" for i in range(ncell)]),
            var=pd.DataFrame(index=genes_all),
        )
        adata.obs.loc["cell0", "splinefit_pseudotime"] = np.nan  # test NaN drop
    except ImportError:
        import pandas as pd

        class _MiniAnnData:  # minimal duck-typed stand-in for the demo
            def __init__(self, X, obs, var_names):
                self.X = X
                self.obs = obs
                self.var_names = var_names
                self.obs_names = obs.index
                self.layers = {}
                self.raw = None
        obs = pd.DataFrame({"splinefit_pseudotime": tt},
                           index=[f"cell{i}" for i in range(ncell)])
        obs.loc["cell0", "splinefit_pseudotime"] = np.nan
        adata = _MiniAnnData(Xcells, obs, genes_all)

    plot_order = ["Lgr5", "Mki67", "Lyz1", "Muc2"]
    res = loess_smoothing_adata(adata, plot_order, span=0.5)
    print(res)
    print("smoothed frame:", res.smoothed.shape,
          "| columns in requested order:", list(res.smoothed.columns))

    fig, _ = plot_loess_curves(res, layout="grid", ncols=2,
                               save="/home/claude/loess_demo_grid.png")
    fig2, _ = plot_loess_curves(res, layout="overlay", zscore=True,
                                save="/home/claude/loess_demo_overlay.png")
    print("saved demo figures.")
