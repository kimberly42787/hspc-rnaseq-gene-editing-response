
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_library_size_qc(
    counts,
    metadata,
    group_col,
    sample_id_col,
    color_map,
    title_prefix="",
    n_sd=3,
    save_path=None,
):
    """
    Plot library size distribution and per-sample library sizes, colored by group.

    Parameters
    ----------
    counts : pd.DataFrame
        Genes x samples counts matrix.
    metadata : pd.DataFrame
        Sample metadata, must contain `group_col` and `sample_id_col`.
    group_col : str
        Column in `metadata` to group/color samples by (e.g. "condition").
    sample_id_col : str
        Column in `metadata` containing sample identifiers matching
        `counts` column names (e.g. "sample_id").
    color_map : dict
        Maps each value in `group_col` to a color, e.g.
        {"RNP_NEG": "#5DCAA5", "GE": "#F2A623", "GE_ANAK": "#E24B4A"}.
    title_prefix : str, optional
        Prefix for plot titles (e.g. "24h").
    n_sd : float, optional
        Number of standard deviations from the mean to flag as an outlier.
    save_path : str or Path, optional
        If given, save the figure to this path.

    Returns
    -------
    library_sizes : pd.Series
        Total counts per sample.
    outliers : pd.Series
        Subset of `library_sizes` flagged as outliers.
    """
    library_sizes = counts.sum(axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].hist(library_sizes / 1e6, bins=15, edgecolor="black")
    axes[0].set_xlabel("Library size (millions of reads)")
    axes[0].set_ylabel("Number of samples")
    axes[0].set_title(f"{title_prefix}: library size distribution".strip(": "))

    sample_order = metadata.sort_values(group_col)[sample_id_col]
    colors = metadata.set_index(sample_id_col).loc[sample_order, group_col].map(color_map)

    axes[1].bar(
        range(len(sample_order)),
        library_sizes.loc[sample_order] / 1e6,
        color=colors,
    )
    axes[1].set_xlabel(f"Sample (grouped by {group_col})")
    axes[1].set_ylabel("Library size (millions of reads)")
    axes[1].set_title(f"{title_prefix}: library size by sample".strip(": "))
    axes[1].set_xticks([])

    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

    mean_lib, std_lib = library_sizes.mean(), library_sizes.std()
    outliers = library_sizes[
        (library_sizes < mean_lib - n_sd * std_lib)
        | (library_sizes > mean_lib + n_sd * std_lib)
    ]
    print(f"Outlier samples (>{n_sd} SD from mean): {len(outliers)}")
    if len(outliers) > 0:
        print(outliers)

    return library_sizes, outliers


def plot_cpm_distribution(counts, min_count, title_prefix="", save_path=None):
    """Plot the distribution of mean log2(CPM+1) per gene, with a vertical
    line marking the CPM-equivalent of a given raw-count filtering threshold."""
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1e6
    mean_log_cpm = np.log2(cpm.mean(axis=1) + 1)

    median_lib_size = library_sizes.median()
    cpm_threshold = min_count / (median_lib_size / 1e6)
    log_threshold = np.log2(cpm_threshold + 1)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(mean_log_cpm, bins=100, edgecolor="none")
    ax.axvline(log_threshold, color="red", linestyle="--",
               label=f"min_count={min_count} threshold (log2CPM={log_threshold:.2f})")
    ax.set_xlabel("Mean log2(CPM + 1) across samples")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"{title_prefix}: gene expression distribution".strip(": "))
    ax.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()


def plot_rescue_distribution(rescue_pct_dict, save_path=None):
    """Plot rescue percentage distributions side by side for multiple groups."""
    n_groups = len(rescue_pct_dict)
    fig, axes = plt.subplots(1, n_groups, figsize=(6 * n_groups, 5), sharey=True)
    if n_groups == 1:
        axes = [axes]

    colors = ["#F2A623", "#E24B4A", "#5DCAA5", "#4C72B0"]

    for ax, (label, values), color in zip(axes, rescue_pct_dict.items(), colors):
        ax.hist(values, bins=20, edgecolor="black", color=color)
        ax.axvline(0, color="black", linewidth=1, linestyle="--")
        ax.axvline(values.median(), color="red", linewidth=1.5, label=f"median = {values.median():.1f}%")
        ax.set_xlabel("Rescue %")
        ax.set_title(f"{label} (n={len(values)} genes)")
        ax.legend()

    axes[0].set_ylabel("Number of genes")
    fig.suptitle("Rescue percentage distribution")
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

def plot_pca(coords, pca, metadata, group_col, sample_id_col, color_map, sample_ids_ordered, title_prefix="", label_points=True, save_path=None):
    """Plot PCA coordinates colored by group, optionally labeled by sample."""
    meta_indexed = metadata.set_index(sample_id_col)
    metadata_ordered = meta_indexed.loc[sample_ids_ordered]

    fig, ax = plt.subplots(figsize=(7, 6))
    for group in metadata_ordered[group_col].unique():
        mask = (metadata_ordered[group_col] == group).values
        ax.scatter(coords[mask, 0], coords[mask, 1], label=group, color=color_map[group], s=80)

    if label_points:
        for i, sample_id in enumerate(sample_ids_ordered):
            label = sample_id.split(":")[1] if ":" in sample_id else sample_id
            ax.annotate(label, (coords[i, 0], coords[i, 1]), fontsize=8, xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title(f"{title_prefix}: PCA by {group_col}".strip(": "))
    ax.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

def plot_volcano(results_df, padj_threshold=0.05, log2fc_threshold=None, title="", save_path=None):
    """Plot a volcano plot (log2FC vs -log10 padj) with significance thresholds marked."""
    df = results_df.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    if log2fc_threshold is not None:
        sig = (df["padj"] < padj_threshold) & (df["log2FoldChange"].abs() > log2fc_threshold)
        sig_label = f"padj<{padj_threshold} & |log2FC|>{log2fc_threshold}"
    else:
        sig = df["padj"] < padj_threshold
        sig_label = f"padj<{padj_threshold}"

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "neg_log10_padj"], s=8, alpha=0.3, color="gray", label="Not significant")
    ax.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "neg_log10_padj"], s=8, alpha=0.6, color="red", label=sig_label)
    ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=0.8)
    if log2fc_threshold is not None:
        ax.axvline(log2fc_threshold, color="black", linestyle="--", linewidth=0.8)
        ax.axvline(-log2fc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

def plot_ma(results_df, padj_threshold=0.05, title="", save_path=None):
    """Plot an MA plot (mean expression vs log2FC) with significant genes highlighted."""
    df = results_df.copy()
    sig = df["padj"] < padj_threshold

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df.loc[~sig, "baseMean"], df.loc[~sig, "log2FoldChange"], s=8, alpha=0.3, color="gray", label="Not significant")
    ax.scatter(df.loc[sig, "baseMean"], df.loc[sig, "log2FoldChange"], s=8, alpha=0.6, color="red", label=f"padj<{padj_threshold}")
    ax.set_xscale("log")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Mean of normalized counts")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

def plot_gsea_dotplot(gsea_df, top_n=15, title="", save_path=None):
    """Plot top GSEA terms by FDR as a horizontal dot plot, colored by NES."""
    df = gsea_df.copy().sort_values("FDR q-val").head(top_n).sort_values("NES")

    fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.35)))
    colors = ["#4C72B0" if nes < 0 else "#C44E52" for nes in df["NES"]]
    ax.barh(df["Term"], df["NES"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("NES")
    ax.set_title(title)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()