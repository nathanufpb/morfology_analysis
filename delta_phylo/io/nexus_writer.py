"""
NexusWriter: export morphological matrices and trees to NEXUS format.

NEXUS is the standard format used by MrBayes, PAUP*, and many other
phylogenetic software packages.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from delta_phylo.matrix.encoding import MatrixEncoder
from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class NexusWriter:
    """Write a morphological matrix to NEXUS format.

    Args:
        matrix: The morphological matrix to export.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix
        self.encoder = MatrixEncoder(matrix)

    def write(self, filepath: str | Path, interleave: bool = False) -> None:
        """Write the NEXUS file.

        Args:
            filepath: Output file path.
            interleave: If True, write interleaved DATA block (not yet
                implemented; always writes sequential format).
        """
        filepath = Path(filepath)
        content = self._build_nexus()
        filepath.write_text(content, encoding="utf-8")
        logger.info("NEXUS file written to %s", filepath)

    def _build_nexus(self) -> str:
        """Build the NEXUS string.

        Returns:
            NEXUS-formatted string.
        """
        taxa_names = self.matrix.get_taxa_names()
        n_taxa = len(taxa_names)
        n_chars = len(self.matrix.get_character_names())

        # Determine number of states for symbols
        max_state = 0
        arr = self.matrix.to_numpy(missing_value=-1)
        if arr.size > 0:
            valid = arr[arr != -1]
            if valid.size > 0:
                max_state = int(valid.max())

        if max_state <= 9:
            symbols = "0123456789"[: max_state + 1]
        else:
            symbols = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[: max_state + 1]

        lines = ["#NEXUS", "", "BEGIN TAXA;"]
        lines.append(f"\tDIMENSIONS NTAX={n_taxa};")
        lines.append("\tTAXLABELS")
        for name in taxa_names:
            # Quote names with spaces
            safe_name = f"'{name}'" if " " in name else name
            lines.append(f"\t\t{safe_name}")
        lines.append("\t;")
        lines.append("END;")
        lines.append("")
        lines.append("BEGIN CHARACTERS;")
        lines.append(f"\tDIMENSIONS NCHAR={n_chars};")
        lines.append(
            f"\tFORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"{symbols}\";"
        )
        lines.append("\tMATRIX")
        max_name_len = max(len(n) for n in taxa_names)
        for name in taxa_names:
            symbol_str = self.encoder.symbol_string(name)
            safe_name = f"'{name}'" if " " in name else name
            padded = safe_name.ljust(max_name_len + 2)
            lines.append(f"\t\t{padded}  {symbol_str}")
        lines.append("\t;")
        lines.append("END;")
        lines.append("")
        # MrBayes block for morphological Mk model
        lines.append("BEGIN MRBAYES;")
        lines.append("\tset autoclose=yes;")
        lines.append(f"\tlset coding=all rates=gamma ngammacat=4;")
        lines.append(f"\tmcmc ngen=1000000 printfreq=1000 samplefreq=500;")
        lines.append("\tsumt burnin=250;")
        lines.append("END;")
        lines.append("")
        return "\n".join(lines)
