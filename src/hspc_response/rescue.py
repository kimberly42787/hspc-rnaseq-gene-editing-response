
import pandas as pd

def compute_rescue_percent(ge_results, ge_anak_results, rnp_neg_baseline=0.0, value_col="log2FoldChange"):

    """
    Compute a continuous rescue percentage per gene.

    rescue % = (GE_value - GE_ANAK_value) / (GE_value - baseline) * 100

    where GE_value is the log2FoldChange of GE vs RNP_NEG (i.e., how far GE
    moved from baseline), and GE_ANAK_value is the log2FoldChange of
    GE_ANAK vs RNP_NEG (how far GE_ANAK moved from the same baseline).

    A value near 100 means GE_ANAK returned close to baseline (full rescue).
    A value near 0 means GE_ANAK is as far from baseline as GE (persistent).
    Negative or >100 values indicate overshoot or reversal, and are kept
    as-is rather than clipped, since they are informative edge cases.

    Parameters
    ----------
    ge_results : pd.DataFrame
        Results table for GE_vs_RNP_NEG, indexed by gene, with `value_col`.
    ge_anak_results : pd.DataFrame
        Results table for GE_ANAK_vs_RNP_NEG, indexed by gene, with `value_col`.
    rnp_neg_baseline : float, optional
        The baseline value (0.0, since log2FC is already relative to RNP_NEG).
    value_col : str, optional
        Column to use for the fold-change values.

    Returns
    -------
    pd.Series
        Rescue percentage per gene, indexed by gene symbol.
    """
    common_genes = ge_results.index.intersection(ge_anak_results.index)
    ge_vals = ge_results.loc[common_genes, value_col]
    ge_anak_vals = ge_anak_results.loc[common_genes, value_col]

    denominator = ge_vals - rnp_neg_baseline
    rescue_pct = (ge_vals - ge_anak_vals) / denominator * 100

    return rescue_pct

def classify_rescue(ge_anak_vs_ge_results, ge_anak_vs_rnp_neg_results, alpha=0.05):
    """
    Classify each gene's rescue status using formal significance tests.

    - "full_rescue": GE_ANAK significantly differs from GE (real movement)
      AND does not significantly differ from RNP_NEG (returned to baseline)
    - "partial_rescue": GE_ANAK significantly differs from GE AND still
      significantly differs from RNP_NEG (moved, but not fully back)
    - "persistent": GE_ANAK does not significantly differ from GE

    Parameters
    ----------
    ge_anak_vs_ge_results : pd.DataFrame
        Results table for GE_ANAK_vs_GE contrast, with a 'padj' column.
    ge_anak_vs_rnp_neg_results : pd.DataFrame
        Results table for GE_ANAK_vs_RNP_NEG contrast, with a 'padj' column.
    alpha : float, optional
        Significance threshold.

    Returns
    -------
    pd.Series
        Categorical rescue classification per gene.
    """
    common_genes = ge_anak_vs_ge_results.index.intersection(ge_anak_vs_rnp_neg_results.index)

    moved_from_ge = ge_anak_vs_ge_results.loc[common_genes, "padj"] < alpha
    still_differs_from_rnp = ge_anak_vs_rnp_neg_results.loc[common_genes, "padj"] < alpha

    classification = pd.Series(index=common_genes, dtype="object")
    classification[moved_from_ge & ~still_differs_from_rnp] = "full_rescue"
    classification[moved_from_ge & still_differs_from_rnp] = "partial_rescue"
    classification[~moved_from_ge] = "persistent"

    return classification

