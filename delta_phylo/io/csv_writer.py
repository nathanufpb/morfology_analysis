"""
CSVWriter: export morphological matrices to CSV format.

CSV export is useful for exploratory data analysis in R, Excel, or other tools.
"""

from __future__ import annotations

import logging
from pathlib import Path

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class CSVWriter:
    """Write a morphological matrix to CSV format.

    Args:
        matrix: The morphological matrix to export.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix

    def write(self, filepath: str | Path, na_rep: str = "?") -> None:
        """Write the CSV file.

        Args:
            filepath: Output file path.
            na_rep: String representation for missing data. Default ``"?"``.
        """
        filepath = Path(filepath)
        self.matrix.df.to_csv(filepath, na_rep=na_rep, index_label="taxon")
        logger.info("CSV file written to %s", filepath)
