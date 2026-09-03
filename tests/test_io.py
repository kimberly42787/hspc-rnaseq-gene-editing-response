"""
Tests for src/hspc_response/io.py
"""

from pathlib import Path
import pytest
import pandas as pd

from hspc_response.io import find_repo_root, load_config, resolve_path, parse_sample_id, build_metadata, load_counts_matrix

def test_find_repo_root_returns_a_directory():
    """The repo root should be a real, existing directory."""
    root = find_repo_root()
    assert root.exists()
    assert root.is_dir()


def test_find_repo_root_contains_expected_project_files():
    """The repo root should contain files we know exist there."""
    root = find_repo_root()
    assert (root / "pyproject.toml").exists()
    assert (root / "config" / "config.yaml").exists()
    assert (root / "src" / "hspc_response").is_dir()


def test_find_repo_root_raises_if_marker_not_found():
    """A marker that doesn't exist anywhere should raise, not fail silently."""
    with pytest.raises(FileNotFoundError):
        find_repo_root(marker="this_file_definitely_does_not_exist.xyz")


def test_load_config_returns_a_dict():
    config = load_config()
    assert isinstance(config, dict)


def test_load_config_has_expected_top_level_keys():
    """Sanity check that the config has the sections we built."""
    config = load_config()
    assert "paths" in config
    assert "design" in config
    assert "thresholds" in config


def test_resolve_path_returns_absolute_path():
    resolved = resolve_path("data/raw")
    assert resolved.is_absolute()


def test_resolve_path_matches_repo_root_plus_relative():
    """resolve_path should equal repo_root / the given relative path."""
    root = find_repo_root()
    resolved = resolve_path("config/config.yaml")
    assert resolved == root / "config" / "config.yaml"


def test_parse_sample_id_handles_numeric_well_code():
    """this verifies if the function parse through a sample ID with numeric-style well code"""
    result = parse_sample_id("BAM:H2:GE_ANAK_24h")
    assert result == {
        "sample_id" : "BAM:H2:GE_ANAK_24h",
        "well_code" : "H2",
        "condition" : "GE_ANAK",
        "timepoint" : "24h"
    }

def test_parse_sample_id_handles_letter_well_code():
    """this verifies if the function parse through a sample ID with letter + numeric-style well code"""
    result = parse_sample_id("BAM:G3:GE_24h")
    assert result == {
        "sample_id" : "BAM:G3:GE_24h",
        "well_code" : "G3",
        "condition" : "GE",
        "timepoint" : "24h"
    }


def test_parse_sample_id_handles_ge_anak_condition():
    """Regression test: GE_ANAK must not be misparsed as GE (this was the
    exact bug risk the regex ordering was designed to avoid)."""
    result = parse_sample_id("BAM:H2:GE_ANAK_24h")
    assert result["condition"] == "GE_ANAK"


def test_parse_sample_id_raises_on_unrecognized_pattern():
    with pytest.raises(ValueError):
        parse_sample_id("this_is_not_a_valid_sample_id")


def test_build_metadata_returns_correct_number_of_rows():
    fake_counts = pd.DataFrame({
        "BAM:N1:GE_96h": [1, 2, 3],
        "BAM:A3:RNP_NEG_24h": [4, 5, 6],
        "BAM:H2:GE_ANAK_24h": [7, 8, 9],
    })
    result = build_metadata(fake_counts)
    assert len(result) == 3


def test_build_metadata_has_expected_columns():
    fake_counts = pd.DataFrame({
        "BAM:N1:GE_96h": [1, 2, 3],
        "BAM:A3:RNP_NEG_24h": [4, 5, 6],
    })
    result = build_metadata(fake_counts)
    assert set(result.columns) == {"sample_id", "well_code", "condition", "timepoint"}


def test_load_counts_matrix_reads_whitespace_separated_file(tmp_path):
    # Create a small fake featureCounts-style file
    fake_file = tmp_path / "fake_counts.txt"
    fake_file.write_text(
        "sample_1 sample_2\n"
        "GENE_A 10 20\n"
        "GENE_B 30 40\n"
    )

    result = load_counts_matrix(fake_file)

    assert result.shape == (2, 2)
    assert list(result.columns) == ["sample_1", "sample_2"]
    assert result.loc["GENE_A", "sample_1"] == 10