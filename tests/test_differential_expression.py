import pandas as pd
import numpy as np
from hspc_response.differential_expression import build_dds, run_contrast, flag_nonconverged


def make_synthetic_data():
    """Small synthetic counts matrix: 20 genes, 6 samples, 2 conditions,
    with a few genes designed to have a clear expression difference."""
    np.random.seed(42)
    n_genes = 20
    samples = ["sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6"]
    conditions = ["control", "control", "control", "treated", "treated", "treated"]

    # Baseline counts around 100, with 5 genes strongly upregulated in "treated"
    counts = np.random.poisson(100, size=(n_genes, 6))
    counts[:5, 3:] = np.random.poisson(500, size=(5, 3))  # genes 0-4 upregulated in treated

    counts_df = pd.DataFrame(counts, columns=samples, index=[f"gene_{i}" for i in range(n_genes)])
    metadata_df = pd.DataFrame({"sample_id": samples, "condition": conditions})

    return counts_df, metadata_df


def test_build_dds_returns_fitted_object():
    counts, metadata = make_synthetic_data()
    dds = build_dds(
        counts, metadata,
        group_col="condition", sample_id_col="sample_id",
        design="~condition", reference_level="control",
    )
    # A fitted dds should have dispersions and LFC convergence info available
    assert "_LFC_converged" in dds.var.columns


def test_run_contrast_detects_designed_effect():
    counts, metadata = make_synthetic_data()
    dds = build_dds(
        counts, metadata,
        group_col="condition", sample_id_col="sample_id",
        design="~condition", reference_level="control",
    )
    result = run_contrast(dds, group_col="condition", numerator="treated", denominator="control")

    # The 5 genes deliberately upregulated in "treated" should show positive log2FC
    upregulated_genes = [f"gene_{i}" for i in range(5)]
    assert (result.loc[upregulated_genes, "log2FoldChange"] > 0).all()


def test_flag_nonconverged_adds_converged_column():
    counts, metadata = make_synthetic_data()
    dds = build_dds(
        counts, metadata,
        group_col="condition", sample_id_col="sample_id",
        design="~condition", reference_level="control",
    )
    result = run_contrast(dds, group_col="condition", numerator="treated", denominator="control")
    flagged = flag_nonconverged(dds, result)

    assert "converged" in flagged.columns