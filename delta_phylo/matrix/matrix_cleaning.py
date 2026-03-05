"""
MatrixCleaner: utilities for cleaning morphological matrices.

Provides methods to remove uninformative characters, impute missing data,
and filter taxa/characters based on completeness thresholds.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


class MatrixCleaner:
    """Clean and filter a :class:`MorphologicalMatrix`.

    Args:
        matrix: The matrix to clean.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix

    def remove_constant_characters(self) -> MorphologicalMatrix:
        """Remove characters that have the same state in all (non-missing) taxa.

        Returns:
            New MorphologicalMatrix without constant characters.
        """
        df = self.matrix.df
        # A column is constant if all non-NaN values are identical
        keep_cols = []
        for col in df.columns:
            non_na = df[col].dropna()
            if non_na.nunique() > 1:
                keep_cols.append(col)

        removed = len(df.columns) - len(keep_cols)
        logger.info("Removing %d constant character(s).", removed)

        new_df = df[keep_cols].copy()
        # Filter character list to match kept columns
        kept_ids = {int(c[1:]) for c in keep_cols}  # strip leading 'C'
        new_chars = [c for c in self.matrix.characters if c.char_id in kept_ids]
        return MorphologicalMatrix(
            df=new_df, characters=new_chars, taxa=self.matrix.taxa
        )

    def remove_autapomorphies(self) -> MorphologicalMatrix:
        """Remove characters that are autapomorphic (unique derived state in one taxon).

        A character is classified as an autapomorphy when **every** derived state
        (state index > 0) it possesses is observed in only a single taxon.  Such
        characters contain no grouping information for cladistic analysis.

        Characters that are constant (same state in all taxa, handled separately
        by :meth:`remove_constant_characters`) are never removed here.

        Returns:
            New MorphologicalMatrix without autapomorphic characters.
        """
        df = self.matrix.df
        keep_cols = []
        for col in df.columns:
            non_na = df[col].dropna()
            if non_na.empty:
                keep_cols.append(col)
                continue
            value_counts = non_na.value_counts()
            # Check whether any derived state (state index > 0) exists at all.
            has_derived = any(state != 0 for state in value_counts.index)
            if not has_derived:
                # Constant at the ancestral state — not autapomorphic; keep it.
                keep_cols.append(col)
                continue
            # Keep if at least one derived state is shared by ≥2 taxa (synapomorphy).
            derived_shared = any(
                int(count) >= 2
                for state, count in value_counts.items()
                if state != 0
            )
            if derived_shared:
                keep_cols.append(col)

        removed = len(df.columns) - len(keep_cols)
        logger.info("Removing %d autapomorphic character(s).", removed)
        new_df = df[keep_cols].copy()
        kept_ids = {int(c[1:]) for c in keep_cols}
        new_chars = [c for c in self.matrix.characters if c.char_id in kept_ids]
        return MorphologicalMatrix(
            df=new_df, characters=new_chars, taxa=self.matrix.taxa
        )

    def filter_by_completeness(
        self,
        min_taxon_completeness: float = 0.5,
        min_char_completeness: float = 0.5,
    ) -> MorphologicalMatrix:
        """Remove taxa and characters below completeness thresholds.

        Args:
            min_taxon_completeness: Minimum fraction of characters scored
                for a taxon to be retained (0–1).
            min_char_completeness: Minimum fraction of taxa scored for a
                character to be retained (0–1).

        Returns:
            Filtered MorphologicalMatrix.
        """
        df = self.matrix.df.copy()
        n_chars = df.shape[1]
        n_taxa = df.shape[0]

        # Filter characters first
        char_completeness = df.notna().mean(axis=0)
        keep_chars = char_completeness[char_completeness >= min_char_completeness].index
        df = df[keep_chars]

        # Filter taxa
        taxon_completeness = df.notna().mean(axis=1)
        keep_taxa = taxon_completeness[taxon_completeness >= min_taxon_completeness].index
        df = df.loc[keep_taxa]

        removed_chars = n_chars - len(keep_chars)
        removed_taxa = n_taxa - len(keep_taxa)
        logger.info(
            "Completeness filter: removed %d characters, %d taxa.",
            removed_chars,
            removed_taxa,
        )

        kept_char_ids = {int(c[1:]) for c in keep_chars}
        kept_taxon_names = set(keep_taxa)
        new_chars = [c for c in self.matrix.characters if c.char_id in kept_char_ids]
        new_taxa = [t for t in self.matrix.taxa if t.name in kept_taxon_names]
        return MorphologicalMatrix(df=df, characters=new_chars, taxa=new_taxa)

    def impute_missing(self, strategy: str = "mode") -> MorphologicalMatrix:
        """Fill missing data using a simple imputation strategy.

        Args:
            strategy: ``"mode"`` fills with the most common state per column.
                      ``"zero"`` fills with 0.

        Returns:
            MorphologicalMatrix with imputed values.
        """
        df = self.matrix.df.copy()
        if strategy == "mode":
            for col in df.columns:
                mode_val = df[col].mode()
                if not mode_val.empty:
                    df[col] = df[col].fillna(mode_val.iloc[0])
        elif strategy == "zero":
            df = df.fillna(0)
        else:
            raise ValueError(f"Unknown imputation strategy: {strategy!r}")

        logger.info("Imputed missing data using strategy=%r.", strategy)
        return MorphologicalMatrix(
            df=df, characters=self.matrix.characters, taxa=self.matrix.taxa
        )
