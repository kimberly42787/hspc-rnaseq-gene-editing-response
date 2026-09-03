

import gseapy as gp
import pandas as pd


def run_gsea_prerank(results_df, gene_sets, ranking_col="log2FoldChange", min_size=15, max_size=500, permutation_num=1000, seed=42):
    """
    Run GSEA prerank on a DESeq2 results table.

    Parameters
    ----------
    results_df : pd.DataFrame
        DESeq2 results, indexed by gene symbol, containing `ranking_col`.
    gene_sets : str or list
        Gene set library name(s) for gseapy (e.g. "MSigDB_Hallmark_2020"),
        or a path to a custom .gmt file.
    ranking_col : str, optional
        Column to rank genes by, default "log2FoldChange" (matches Di Micco
        et al.'s method: "all genes were ranked by decreasing Log2FC values").
    min_size, max_size : int, optional
        Gene set size bounds.
    permutation_num : int, optional
        Number of permutations for significance testing.
    seed : int, optional
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        GSEA results table.
    """
    ranked = results_df[[ranking_col]].dropna()
    ranked = ranked.sort_values(ranking_col, ascending=False)
    ranked_series = ranked[ranking_col]

    pre_res = gp.prerank(
        rnk=ranked_series,
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        outdir=None,
        seed=seed,
    )
    return pre_res.res2d