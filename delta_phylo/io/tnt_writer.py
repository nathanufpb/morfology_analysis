"""
TNTWriter: export morphological matrices to TNT (Tree analysis using New Technology) format.

TNT is a fast parsimony program by Goloboff et al. (https://www.lillo.org.ar/phylogeny/tnt/).
"""

from __future__ import annotations

import logging
from pathlib import Path

from delta_phylo.matrix.encoding import MatrixEncoder
from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class TNTWriter:
    """Write a morphological matrix to TNT format.

    Args:
        matrix: The morphological matrix to export.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix
        self.encoder = MatrixEncoder(matrix)

    def write(self, filepath: str | Path) -> None:
        """Write the TNT file.

        Args:
            filepath: Output file path.
        """
        filepath = Path(filepath)
        content = self._build_tnt()
        filepath.write_text(content, encoding="utf-8")
        logger.info("TNT file written to %s", filepath)

    def _build_tnt(self) -> str:
        """Build the TNT string.

        TNT format:
          ``xread``
          ``'<title>'``
          ``<n_chars> <n_taxa>``
          ``<taxon_name>  <sequence>``
          ``;``
          ``proc /;``

        Returns:
            TNT-formatted string.
        """
        taxa_names = self.matrix.get_taxa_names()
        n_taxa = len(taxa_names)
        n_chars = len(self.matrix.get_character_names())

        lines = [
            "xread",
            f"'Morphological matrix exported by delta-phylo'",
            f"{n_chars} {n_taxa}",
        ]

        for name in taxa_names:
            # TNT names cannot contain spaces; replace with underscore
            safe_name = name.replace(" ", "_")
            symbol_str = self.encoder.symbol_string(name)
            lines.append(f"{safe_name}  {symbol_str}")

        lines.append(";")
        lines.append("proc /;")
        lines.append("")
        return "\n".join(lines)
