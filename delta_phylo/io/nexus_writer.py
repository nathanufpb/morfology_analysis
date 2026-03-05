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

    def write(self, filepath: str | Path, interleave: bool = False, block_width: int = 60) -> None:
        """Write the NEXUS file.

        Args:
            filepath: Output file path.
            interleave: If True, write the MATRIX in interleaved blocks.
                Each block contains *block_width* characters per line.
            block_width: Number of characters per taxon per interleaved block
                (default 60). Ignored when *interleave* is False.
        """
        filepath = Path(filepath)
        content = self._build_nexus(interleave=interleave, block_width=block_width)
        filepath.write_text(content, encoding="utf-8")
        logger.info("NEXUS file written to %s", filepath)

    def _build_nexus(self, interleave: bool = False, block_width: int = 60) -> str:
        """Build the NEXUS string.

        Args:
            interleave: Whether to write the MATRIX in interleaved format.
            block_width: Characters per taxon per interleaved block.

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
        interleave_flag = " INTERLEAVE=YES" if interleave else ""
        lines.append(
            f"\tFORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"{symbols}\"{interleave_flag};"
        )
        lines.append("\tMATRIX")
        max_name_len = max(len(n) for n in taxa_names)
        if interleave:
            # Write MATRIX in blocks of block_width characters
            for block_start in range(0, n_chars, block_width):
                block_end = min(block_start + block_width, n_chars)
                lines.append(f"\t\t[ characters {block_start + 1}-{block_end} ]")
                for name in taxa_names:
                    symbol_str = self.encoder.symbol_string(name)
                    block_str = symbol_str[block_start:block_end]
                    safe_name = f"'{name}'" if " " in name else name
                    padded = safe_name.ljust(max_name_len + 2)
                    lines.append(f"\t\t{padded}  {block_str}")
                lines.append("")   # blank line between blocks
        else:
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
