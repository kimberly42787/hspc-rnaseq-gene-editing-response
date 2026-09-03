import pandas as pd
from hspc_response.plotting import plot_library_size_qc


def test_plot_library_size_qc_returns_correct_library_sizes():
    counts = pd.DataFrame({
        "sample_1": [10, 20, 30],
        "sample_2": [5, 15, 25],
    })
    metadata = pd.DataFrame({
        "sample_id": ["sample_1", "sample_2"],
        "condition": ["GE", "RNP_NEG"],
    })
    color_map = {"GE": "#000000", "RNP_NEG": "#FFFFFF"}

    library_sizes, outliers = plot_library_size_qc(
        counts, metadata,
        group_col="condition", sample_id_col="sample_id",
        color_map=color_map,
    )

    assert library_sizes["sample_1"] == 60  # 10+20+30
    assert library_sizes["sample_2"] == 45  # 5+15+25


def test_plot_library_size_qc_flags_known_outlier():
    counts = pd.DataFrame({
        "sample_1": [100, 100, 100],
        "sample_2": [100, 100, 100],
        "sample_3": [100, 100, 100],
        "sample_4": [100, 100, 100],
        "sample_5": [10000, 10000, 10000],  # obvious outlier
    })
    metadata = pd.DataFrame({
        "sample_id": ["sample_1", "sample_2", "sample_3", "sample_4", "sample_5"],
        "condition": ["GE", "GE", "RNP_NEG", "RNP_NEG", "GE_ANAK"],
    })
    color_map = {"GE": "#000000", "RNP_NEG": "#FFFFFF", "GE_ANAK": "#FF0000"}

    library_sizes, outliers = plot_library_size_qc(
        counts, metadata,
        group_col="condition", sample_id_col="sample_id",
        color_map=color_map, n_sd=1,
    )

    assert "sample_5" in outliers.index