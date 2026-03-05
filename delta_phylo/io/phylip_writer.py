"""
PhylipWriter: export morphological matrices to PHYLIP format.

The PHYLIP format is used by the PHYLIP package (Felsenstein 1989+)
and is accepted by many other phylogenetic programs.
"""

from __future__ import annotations

import logging
from pathlib import Path

from delta_phylo.matrix.encoding import MatrixEncoder
from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class PhylipWriter:
    """Write a morphological matrix to PHYLIP format.

    Args:
        matrix: The morphological matrix to export.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix
        self.encoder = MatrixEncoder(matrix)

    def write(self, filepath: str | Path) -> None:
        """Write the PHYLIP file.

        Args:
            filepath: Output file path.
        """
        filepath = Path(filepath)
        content = self._build_phylip()
        filepath.write_text(content, encoding="utf-8")
        logger.info("PHYLIP file written to %s", filepath)

    def _build_phylip(self) -> str:
        """Build the PHYLIP string.

        PHYLIP sequential format:
          Line 1: ``<n_taxa> <n_chars>``
          Subsequent lines: ``<name (10 chars)><sequence>``

        Returns:
            PHYLIP-formatted string.
        """
        taxa_names = self.matrix.get_taxa_names()
        n_taxa = len(taxa_names)
        n_chars = len(self.matrix.get_character_names())

        lines = [f" {n_taxa} {n_chars}"]
        for name in taxa_names:
            # PHYLIP name field is exactly 10 characters (padded or truncated)
            padded_name = name[:10].ljust(10)
            symbol_str = self.encoder.symbol_string(name)
            lines.append(f"{padded_name}{symbol_str}")

        lines.append("")
        return "\n".join(lines)
