"""
MorphologicalMatrix and MatrixBuilder: construct a taxa × characters matrix
from parsed DELTA objects.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

from delta_phylo.parser.characters import Character
from delta_phylo.parser.taxa import Taxon

logger = logging.getLogger(__name__)

# Sentinel value for missing data in numeric arrays
MISSING_INT: int = -1


@dataclass
class MorphologicalMatrix:
    """Stores a morphological data matrix as a Pandas DataFrame.

    Rows correspond to taxa, columns to characters.
    Missing data is represented as ``NaN`` (float) or ``-1`` (int arrays).

    Attributes:
        df: The underlying DataFrame with taxa names as index and character
            IDs/names as columns.
        characters: Ordered list of Character objects matching columns.
        taxa: Ordered list of Taxon objects matching rows.
    """

    df: pd.DataFrame
    characters: List[Character] = field(default_factory=list)
    taxa: List[Taxon] = field(default_factory=list)

    def __repr__(self) -> str:
        return (
            f"MorphologicalMatrix(taxa={len(self.taxa)}, "
            f"characters={len(self.characters)})"
        )

    @property
    def shape(self):
        """Return (num_taxa, num_characters) shape tuple."""
        return self.df.shape

    def to_numpy(self, missing_value: int = MISSING_INT) -> np.ndarray:
        """Return the matrix as a NumPy integer array.

        Polymorphic states are collapsed to the *first* state.
        **Note**: numeric character values are truncated to int; use
        :meth:`to_numpy_float` to preserve them.

        Args:
            missing_value: Integer sentinel for missing data. Default ``-1``.

        Returns:
            2D integer array of shape (n_taxa, n_chars).
        """
        arr = self.df.values.copy()
        # Replace NaN with missing_value
        result = np.where(np.isnan(arr.astype(float)), missing_value, arr)
        return result.astype(int)

    def to_numpy_float(self, missing_value: float = float("nan")) -> np.ndarray:
        """Return the matrix as a NumPy float64 array.

        Unlike :meth:`to_numpy`, this method preserves decimal values for
        numeric (continuous) characters.

        Args:
            missing_value: Sentinel for missing data (default ``nan``).

        Returns:
            2D float64 array of shape (n_taxa, n_chars).
        """
        arr = self.df.values.copy().astype(float)
        if not np.isnan(missing_value):
            arr = np.where(np.isnan(arr), missing_value, arr)
        return arr

    def get_taxa_names(self) -> List[str]:
        """Return list of taxon names (row index)."""
        return list(self.df.index)

    def get_character_names(self) -> List[str]:
        """Return list of character column names."""
        return list(self.df.columns)


class MatrixBuilder:
    """Build a :class:`MorphologicalMatrix` from parsed DELTA objects.

    Usage::

        builder = MatrixBuilder(characters, taxa)
        matrix = builder.build()

    Args:
        characters: Ordered list of Character objects.
        taxa: Ordered list of Taxon objects.
        polymorphic_strategy: How to handle polymorphic states.
            - ``"first"``    – use the first (lowest/minimum) state.
            - ``"missing"``  – treat polymorphic as missing data.
            - ``"majority"`` – use the most common state (mode); ties are
              broken by taking the smallest tied state.
    """

    _POLYMORPHIC_STRATEGIES = ("first", "missing", "majority")

    def __init__(
        self,
        characters: List[Character],
        taxa: List[Taxon],
        polymorphic_strategy: str = "first",
    ) -> None:
        if polymorphic_strategy not in self._POLYMORPHIC_STRATEGIES:
            raise ValueError(
                f"polymorphic_strategy must be one of {self._POLYMORPHIC_STRATEGIES}"
            )
        self.characters = characters
        self.taxa = taxa
        self.polymorphic_strategy = polymorphic_strategy

    def build(self) -> MorphologicalMatrix:
        """Construct and return the morphological matrix.

        Returns:
            A populated :class:`MorphologicalMatrix` instance.
        """
        logger.info(
            "Building matrix: %d taxa × %d characters",
            len(self.taxa),
            len(self.characters),
        )
        col_names = [f"C{c.char_id}" for c in self.characters]
        char_ids = [c.char_id for c in self.characters]
        rows: Dict[str, List[Optional[float]]] = {}

        for taxon in self.taxa:
            row: List[Optional[float]] = []
            for char_id in char_ids:
                score = taxon.get_score(char_id)
                row.append(self._resolve_score(score))
            rows[taxon.name] = row

        df = pd.DataFrame.from_dict(rows, orient="index", columns=col_names)
        df = df.astype(float)
        return MorphologicalMatrix(df=df, characters=self.characters, taxa=self.taxa)

    def _resolve_score(self, score) -> Optional[float]:
        """Convert a raw score to a float suitable for the DataFrame.

        For numeric (continuous) characters the float value is stored as-is.
        For discrete characters the value is converted to a 0-based state index.

        Args:
            score: int, float, list[int], or None.

        Returns:
            Float value or NaN for missing data.
        """
        if score is None:
            return float("nan")
        # Float values come from NUMERIC characters — store directly.
        if isinstance(score, float):
            return score
        if isinstance(score, list):
            if not score:
                return float("nan")
            if self.polymorphic_strategy == "first":
                return float(min(score))
            if self.polymorphic_strategy == "majority":
                # Mode: most common state; ties broken by smallest state.
                counts: Dict[int, int] = {}
                for s in score:
                    counts[s] = counts.get(s, 0) + 1
                max_count = max(counts.values())
                modal_states = sorted(k for k, v in counts.items() if v == max_count)
                return float(modal_states[0])
            # "missing" strategy
            return float("nan")
        return float(score)
