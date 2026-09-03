import pandas as pd
from hspc_response.qc import filter_low_gene_expression


def test_filter_low_gene_expression_removes_low_count_genes():
    counts = pd.DataFrame({
        "sample_1": [1000, 0, 1000],
        "sample_2": [1000, 0, 1000],
        "sample_3": [1000, 0, 1000],
    }, index=["HIGH_GENE", "LOW_GENE", "HIGH_GENE_2"])

    metadata = pd.DataFrame({
        "sample_id": ["sample_1", "sample_2", "sample_3"],
        "condition": ["GE", "RNP_NEG", "GE_ANAK"],
    })

    counts_filtered, summary = filter_low_gene_expression(counts, metadata, group_col="condition")

    assert "HIGH_GENE" in counts_filtered.index
    assert "HIGH_GENE_2" in counts_filtered.index
    assert "LOW_GENE" not in counts_filtered.index


def test_filter_low_gene_expression_returns_correct_summary():
    counts = pd.DataFrame({
        "sample_1": [1000, 0],
        "sample_2": [1000, 0],
        "sample_3": [1000, 0],
    }, index=["HIGH_GENE", "LOW_GENE"])

    metadata = pd.DataFrame({
        "sample_id": ["sample_1", "sample_2", "sample_3"],
        "condition": ["GE", "RNP_NEG", "GE_ANAK"],
    })

    counts_filtered, summary = filter_low_gene_expression(counts, metadata, group_col="condition")

    assert summary["genes_before_filtering"] == 2
    assert summary["genes_after_filtering"] == 1


def test_filter_low_gene_expression_respects_smallest_group_size():
    # 6 samples, 2 conditions x 3 replicates each (smallest_group_size = 3).
    # HIGH_IN_HALF passes the CPM threshold in only 3 of 6 samples (all
    # within one condition) -- actually meets smallest_group_size=3, so
    # let's make it fail in enough samples that it CAN'T meet the bar.
    counts = pd.DataFrame({
        "sample_1": [1000, 1000],
        "sample_2": [1000, 0],
        "sample_3": [1000, 0],
        "sample_4": [1000, 0],
        "sample_5": [1000, 0],
        "sample_6": [1000, 0],
    }, index=["ALWAYS_HIGH", "HIGH_IN_ONE"])

    metadata = pd.DataFrame({
        "sample_id": ["sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6"],
        "condition": ["GE", "GE", "GE", "RNP_NEG", "RNP_NEG", "RNP_NEG"],
    })

    counts_filtered, summary = filter_low_gene_expression(counts, metadata, group_col="condition")

    assert "ALWAYS_HIGH" in counts_filtered.index
    assert "HIGH_IN_ONE" not in counts_filtered.index

