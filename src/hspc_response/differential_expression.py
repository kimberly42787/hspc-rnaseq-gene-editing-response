import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA


def build_dds(counts, metadata, group_col, sample_id_col, design, reference_level, n_cpus=4):
    """
    Build a fitted DeseqDataSet for one timepoint's data.

    Parameters
    ----------
    counts : pd.DataFrame
        Genes x samples counts matrix (filtered).
    metadata : pd.DataFrame
        Sample metadata; must contain `group_col` and `sample_id_col`.
    group_col : str
        Column to use as the design factor (e.g. "condition").
    sample_id_col : str
        Column with sample identifiers matching `counts` columns.
    design : str
        DESeq2 design formula (e.g. "~condition").
    reference_level : str
        Reference level for `group_col` (e.g. "RNP_NEG").
    n_cpus : int, optional
        Number of CPUs for DESeq2 inference.

    Returns
    -------
    dds : DeseqDataSet
        Fitted DESeq2 dataset (after calling .deseq2()).
    """
    counts_t = counts.T
    meta_indexed = metadata.set_index(sample_id_col)
    metadata_ordered = meta_indexed.loc[counts_t.index].copy()

    categories = [reference_level] + [
        lvl for lvl in metadata_ordered[group_col].unique() if lvl != reference_level
    ]
    metadata_ordered[group_col] = pd.Categorical(metadata_ordered[group_col], categories=categories)

    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(
        counts=counts_t,
        metadata=metadata_ordered,
        design=design,
        inference=inference,
    )
    dds.deseq2()
    return dds

def run_contrast(dds, group_col, numerator, denominator):
    """
    Extract a single pairwise contrast from a fitted DeseqDataSet.

    Parameters
    ----------
    dds : DeseqDataSet
        A fitted (post-.deseq2()) dataset.
    group_col : str
        The design factor column (e.g. "condition").
    numerator : str
        The level being compared (e.g. "GE").
    denominator : str
        The reference/baseline level (e.g. "RNP_NEG").

    Returns
    -------
    results_df : pd.DataFrame
        DESeq2 results table for this contrast.
    """
    stat = DeseqStats(dds, contrast=[group_col, numerator, denominator])
    stat.summary()
    return stat.results_df

def flag_nonconverged(dds, results_df):
    """
    Add a 'converged' boolean column to a results table, based on the
    fitted dds's LFC convergence flags.
    """
    non_converged_genes = dds.var.index[~dds.var["_LFC_converged"]]
    results_df = results_df.copy()
    results_df["converged"] = ~results_df.index.isin(non_converged_genes)
    return results_df


def compute_pca(dds, n_components=2):
    """Compute PCA on VST-transformed counts from a fitted DeseqDataSet."""
    dds.vst_fit()
    vst_counts = dds.vst_transform()
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(vst_counts)
    return coords, pca