"""
Validation utilities for delta-phylo inputs and intermediate objects.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


def validate_delta_directory(path: str | Path) -> Path:
    """Check that a path points to a readable DELTA dataset directory.

    Args:
        path: Path to the DELTA dataset directory.

    Returns:
        Resolved :class:`pathlib.Path` object.

    Raises:
        FileNotFoundError: If the path does not exist.
        NotADirectoryError: If the path is not a directory.
        ValueError: If no recognised DELTA files are found.
    """
    p = Path(path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"DELTA dataset path not found: {p}")
    if not p.is_dir():
        raise NotADirectoryError(f"Expected a directory, got: {p}")

    recognised = {"items", "characters", "chars", "specs"}
    found = {
        f.name.lower()
        for f in p.iterdir()
        if f.is_file() and f.name.lower() in recognised
    }
    if not found:
        raise ValueError(
            f"No recognised DELTA files (items, characters, specs) found in {p}. "
            "Ensure the directory contains at least an 'items' file."
        )
    logger.debug("Validated DELTA directory: %s (files: %s)", p, found)
    return p


def validate_matrix(matrix: MorphologicalMatrix, min_taxa: int = 2) -> None:
    """Validate that a morphological matrix is suitable for analysis.

    Args:
        matrix: MorphologicalMatrix to validate.
        min_taxa: Minimum number of taxa required.

    Raises:
        ValueError: If validation fails.
    """
    n_taxa, n_chars = matrix.shape
    if n_taxa < min_taxa:
        raise ValueError(
            f"Matrix has only {n_taxa} taxa; at least {min_taxa} required."
        )
    if n_chars == 0:
        raise ValueError("Matrix has no characters.")

    # Check proportion of missing data
    total_cells = matrix.df.size
    missing_cells = matrix.df.isna().sum().sum()
    missing_frac = missing_cells / total_cells if total_cells > 0 else 0.0

    if missing_frac > 0.9:
        logger.warning(
            "Matrix has >90%% missing data (%.1f%%). Results may be unreliable.",
            missing_frac * 100,
        )
    logger.debug(
        "Matrix validated: %d taxa × %d characters (missing: %.1f%%)",
        n_taxa,
        n_chars,
        missing_frac * 100,
    )
