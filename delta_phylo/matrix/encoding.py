"""
MatrixEncoder: encode morphological matrices for various phylogenetic formats.

Handles binary encoding, ordered multistate characters, and gap coding.
"""

from __future__ import annotations

import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class MatrixEncoder:
    """Encode a morphological matrix for phylogenetic analysis.

    Provides methods to convert multistate characters to binary (additive/
    non-additive) and to produce NEXUS-compatible symbol mappings.

    Args:
        matrix: The matrix to encode.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix

    def to_binary(self) -> MorphologicalMatrix:
        """Convert all multistate characters to binary presence/absence columns.

        For a character with *k* states, *k* binary columns are created,
        one per state, following the standard additive binary coding.

        Returns:
            New MorphologicalMatrix with all-binary columns.
        """
        df = self.matrix.df
        new_cols: dict = {}
        new_chars: list = []

        for char in self.matrix.characters:
            col_name = f"C{char.char_id}"
            if col_name not in df.columns:
                continue
            col = df[col_name]
            max_state = int(col.dropna().max()) if not col.dropna().empty else 0

            if max_state <= 1:
                # Already binary
                new_cols[col_name] = col
                new_chars.append(char)
            else:
                # One column per state
                for state_idx in range(max_state + 1):
                    new_col_name = f"C{char.char_id}_s{state_idx}"
                    new_cols[new_col_name] = col.apply(
                        lambda x, s=state_idx: (
                            float("nan") if pd.isna(x) else float(x == s)
                        )
                    )
                    from delta_phylo.parser.characters import Character, CharacterType

                    sub_char = Character(
                        char_id=char.char_id,
                        name=f"{char.name} [state {state_idx}]",
                        char_type=CharacterType.BINARY,
                    )
                    new_chars.append(sub_char)

        new_df = pd.DataFrame(new_cols, index=df.index)
        return MorphologicalMatrix(
            df=new_df, characters=new_chars, taxa=self.matrix.taxa
        )

    def symbol_string(self, taxon_name: str) -> str:
        """Return a NEXUS/TNT state symbol string for a single taxon.

        States are encoded as integers 0–9 or letters A–Z for higher states.
        Missing data is ``?``.

        Args:
            taxon_name: Name of the taxon row.

        Returns:
            String of state symbols (one per character).
        """
        if taxon_name not in self.matrix.df.index:
            raise KeyError(f"Taxon {taxon_name!r} not found in matrix.")

        row = self.matrix.df.loc[taxon_name]
        symbols = []
        for val in row:
            if pd.isna(val):
                symbols.append("?")
            else:
                symbols.append(self._encode_state(int(val)))
        return "".join(symbols)

    @staticmethod
    def _encode_state(state: int) -> str:
        """Encode a state integer as a single character symbol.

        States 0–9 map to ``'0'``–``'9'``, states 10–35 map to ``'A'``–``'Z'``.

        Args:
            state: 0-based integer state.

        Returns:
            Single character symbol.
        """
        if state < 0:
            return "?"
        if state <= 9:
            return str(state)
        letter_idx = state - 10
        if letter_idx < 26:
            return chr(ord("A") + letter_idx)
        # Fallback for very large state spaces
        return "?"
