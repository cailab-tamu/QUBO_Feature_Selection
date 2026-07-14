"""
feature_report.py

Bridge between the QUBO/annealing output and gene names: builds the feature
DataFrame that annealing_functions.quantum_test expects, and reports / plots the
selected features.

Typical flow
------------
    adata_f = prefilter_genes(adata, min_cells_frac=0.01, top_n_variable=5000)
    Xy      = mat_transform_sc(adata_f, obs_key="dpt_pseudotime_rank")
    MI_mat  = mutual_information_matrix(Xy, bins_method="rice", binning="quantile")

    df_features = make_feature_df(adata_f, MI_mat)      # names MUST come from adata_f
    alpha, xsol = annealing_process(MI_mat, K=100)
    new_df, qa_sol = quantum_test(alpha, MI_mat, df_features, K=100, mode="sa")

    print_selected_features(new_df, top_n=30)
    plot_selected_features(adata_f, new_df, pseudotime_key="dpt_pseudotime_rank")
"""

import numpy as np
import pandas as pd


def make_feature_df(adata, MI_mat=None):
    """
    Build the feature DataFrame that quantum_test expects.

    quantum_test does `df_features.rename(columns={0: 'Features'})`, so the gene
    names must live in a column literally named 0 (the default when you do
    pd.DataFrame(list_of_names)).

    CRITICAL - alignment: the gene names must come from the SAME AnnData that
    produced MI_mat. If you pre-filtered genes, pass the FILTERED object
    (adata_f), not the original. Passing MI_mat lets this function verify the
    counts match and fail loudly instead of silently mislabelling every gene.

    Parameters
    ----------
    adata : AnnData
        The object used to build Xy (i.e. post-prefilter).
    MI_mat : ndarray, optional
        The MI matrix. If given, checks that MI_mat.shape[0] == n_genes + 1
        (the +1 is the appended target row).

    Returns
    -------
    df_features : pd.DataFrame, one column named 0 holding gene names, in the
        same order as the rows of MI_mat (excluding the final target row).
    """
    genes = list(adata.var_names)
    df = pd.DataFrame(genes)          # column name is the integer 0

    if MI_mat is not None:
        expected = MI_mat.shape[0] - 1     # last row/col is the target
        if len(genes) != expected:
            raise ValueError(
                f"Gene/MI mismatch: adata has {len(genes)} genes but MI_mat "
                f"implies {expected} features (MI_mat is {MI_mat.shape}).\n"
                f"Did you build MI_mat from a pre-filtered AnnData and then "
                f"pass the UNFILTERED one here? Pass the same object you gave "
                f"to mat_transform_sc."
            )
        print(f"[make_feature_df] {len(genes):,} genes, aligned with MI_mat "
              f"{MI_mat.shape} (last row = target). OK.")
    else:
        print(f"[make_feature_df] {len(genes):,} genes (no MI_mat given, "
              f"alignment NOT verified).")
    return df


def sample_to_vector(sampleset, n_features):
    """
    Extract a selection vector from a dimod SampleSet in VARIABLE-INDEX order.

    annealing_functions.quantum_test uses `list(sample.values())`, which returns
    values in dict insertion order rather than variable order. dimod normally
    yields sorted variables so it usually agrees, but it is not guaranteed; if it
    ever differs, gene i silently receives another gene's bit. This does the
    lookup explicitly.

    Returns
    -------
    ndarray of int (0/1), length n_features.
    """
    sample = sampleset.first.sample
    missing = [i for i in range(n_features) if i not in sample]
    if missing:
        raise KeyError(f"{len(missing)} variable(s) absent from the sample "
                       f"(e.g. {missing[:5]}). n_features={n_features}?")
    return np.array([int(sample[i]) for i in range(n_features)], dtype=int)


def selected_genes(new_df, gene_col="Features"):
    """Return the list of selected gene names from quantum_test's output."""
    if gene_col not in new_df.columns:
        raise KeyError(f"'{gene_col}' not in {list(new_df.columns)}")
    return new_df[gene_col].tolist()


def print_selected_features(new_df, top_n=None, gene_col="Features",
                            score_col="feature_score"):
    """
    Print the QUBO-selected features.

    quantum_test sorts ascending by feature_score, and because the QUBO is a
    MINIMIZATION, more negative = more strongly selected. So the top of the list
    is the most important, do not re-sort descending expecting "best first".

    Parameters
    ----------
    new_df : pd.DataFrame
        First return value of quantum_test.
    top_n : int or None
        Print only the first N (most negative score). None prints all.
    """
    n = len(new_df)
    print(f"QUBO selected {n:,} feature(s)")
    if score_col in new_df.columns:
        s = new_df[score_col]
        print(f"feature_score: min={s.min():.4f}  max={s.max():.4f}  "
              f"(more negative = more strongly selected)")
    print("-" * 46)

    show = new_df if top_n is None else new_df.head(top_n)
    for rank, (_, row) in enumerate(show.iterrows(), start=1):
        name = row[gene_col]
        if score_col in new_df.columns:
            print(f"{rank:>4}. {name:<20} {row[score_col]:>10.4f}")
        else:
            print(f"{rank:>4}. {name}")
    if top_n is not None and n > top_n:
        print(f"... and {n - top_n:,} more")
    return selected_genes(new_df, gene_col)


def save_selected_features(new_df, path="filt_df_QA.csv"):
    """Write the selected-feature table to CSV (matches the README's output)."""
    new_df.to_csv(path, index=False)
    print(f"[save_selected_features] wrote {len(new_df):,} rows -> {path}")
    return path


def plot_selected_features(adata, new_df, pseudotime_key, top_n=9,
                           gene_col="Features", layer=None, span=0.75,
                           ncols=3, save=None):
    """
    LOESS-smooth the top selected genes along pseudotime and plot them.

    Requires loess_smoothing.py to be importable (put it in the qfeatures
    package alongside this module).

    Parameters
    ----------
    adata : AnnData
        Same (filtered) object used for the QUBO, so the gene names resolve.
    new_df : pd.DataFrame
        quantum_test output.
    pseudotime_key : str
        obs column, e.g. "dpt_pseudotime_rank".
    top_n : int
        How many of the top-selected genes to plot.
    """
    try:
        from .loess_smoothing import loess_smoothing_adata, plot_loess_curves
    except ImportError:
        from loess_smoothing import loess_smoothing_adata, plot_loess_curves

    genes = selected_genes(new_df, gene_col)[:top_n]
    genes = [g for g in genes if g in set(adata.var_names)]
    if not genes:
        raise ValueError("None of the selected genes are in adata.var_names. "
                         "Are you passing the filtered AnnData?")

    print(f"[plot_selected_features] plotting {len(genes)} gene(s) vs "
          f"'{pseudotime_key}': {genes}")
    res = loess_smoothing_adata(adata, genes, pseudotime_key=pseudotime_key,
                                layer=layer, span=span)
    return plot_loess_curves(res, layout="grid", ncols=ncols, save=save)
