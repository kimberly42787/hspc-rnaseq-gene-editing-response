
import pandas as pd


def filter_low_gene_expression(counts: pd.DataFrame, metadata: pd.DataFrame, group_col: str, min_count: int = 10) -> tuple[pd.DataFrame, dict]:
    """
    Filter low-expressing genes using a CPM-based threshold.

    A gene is kept if its CPM meets the threshold (equivalent to `min_count`
    raw counts at the median library size) in at least as many samples as
    the smallest group size, based on `group_col` in `metadata`.

    Parameters
    ----------
    counts : pd.DataFrame
        Genes x samples counts matrix.
    metadata : pd.DataFrame
        Sample metadata; must contain `group_col`, with samples matching
        `counts` columns.
    group_col : str
        Column in `metadata` used to determine the smallest group size
        (e.g. "condition").
    min_count : int, optional
        Minimum raw count (at median library size) a gene must reach,
        default 10.

    Returns
    -------
    counts_filtered : pd.DataFrame
        Filtered counts matrix.
    summary : dict
        Before/after gene counts and filtering parameters used.
    """
    library_sizes = counts.sum(axis=0)
    median_lib_size = library_sizes.median()
    cpm_threshold = min_count / (median_lib_size / 1e6)
    cpm = counts.div(library_sizes, axis=1) * 1e6

    smallest_group_size = metadata[group_col].value_counts().min()

    genes_to_keep = (cpm >= cpm_threshold).sum(axis=1) >= smallest_group_size
    counts_filtered = counts.loc[genes_to_keep]

    summary = {
        "genes_before_filtering": counts.shape[0],
        "genes_after_filtering": counts_filtered.shape[0],
        "min_count": min_count,
        "smallest_group_size": smallest_group_size,
    }

    return counts_filtered, summary