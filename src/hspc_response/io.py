"""
Handles loading configuration and locating project files, anchored to the
repo root regardless of the working directory the calling code runs from.
"""

from pathlib import Path
import yaml
import re
import pandas as pd


def find_repo_root(marker: str = "pyproject.toml") -> Path:
    """
    Walk up the directory tree from this file's location until a folder
    containing `marker` is found. Returns that folder as the repo root.

    Raises FileNotFoundError if no such folder exists up to the filesystem root.
    """
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / marker).exists():
            return parent
    raise FileNotFoundError(
        f"Could not find repo root (no '{marker}' found in any parent directory of {current})"
    )


def load_config(config_path: str = "config/config.yaml") -> dict:
    """
    Load the project YAML configuration file.

    `config_path` is interpreted relative to the repo root, not the
    current working directory.
    """
    repo_root = find_repo_root()
    full_path = repo_root / config_path
    with open(full_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def resolve_path(relative_path: str) -> Path:
    """
    Resolve a repo-root-relative path string (e.g. one pulled from the
    config file) into an absolute Path.
    """
    return find_repo_root() / relative_path


def load_counts_matrix(path) -> pd.DataFrame:
    """Load a single featureCounts-style matrix (whitespace-separated,
    gene symbols as index)."""
    return pd.read_csv(path, sep=r"\s+")


def parse_sample_id(sample_id: str) -> dict:
    """Extract well code, condition, and timepoint from a column name
    like 'BAM:N1:GE_96h' or 'BAM:A3:RNP_NEG_24h'."""
    match = re.match(r"BAM:([A-Z]\d+):(RNP_NEG|GE_ANAK|GE)_(\d+h)", sample_id)
    if not match:
        raise ValueError(f"Sample ID did not match expected pattern: {sample_id}")
    well_code, condition, timepoint = match.groups()
    return {
        "sample_id": sample_id,
        "well_code": well_code,
        "condition": condition,
        "timepoint": timepoint,
    }


def build_metadata(counts: pd.DataFrame) -> pd.DataFrame:
    """Build a metadata table by parsing every column name in a counts matrix."""
    rows = [parse_sample_id(col) for col in counts.columns]
    return pd.DataFrame(rows)

