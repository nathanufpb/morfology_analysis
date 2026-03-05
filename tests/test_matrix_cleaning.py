"""Tests for MatrixCleaner — closes issue #7.

Covers:
  - remove_constant_characters()
  - remove_autapomorphies()
  - filter_by_completeness()
  - impute_missing()
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix, MatrixBuilder
from delta_phylo.matrix.matrix_cleaning import MatrixCleaner
from delta_phylo.parser.delta_reader import DeltaReader


# ---------------------------------------------------------------------------
# Helper factories
# ---------------------------------------------------------------------------


def _make_matrix(df: pd.DataFrame) -> MorphologicalMatrix:
    """Wrap a raw DataFrame in a MorphologicalMatrix (no characters/taxa lists)."""
    return MorphologicalMatrix(df=df, characters=[], taxa=[])


def _matrix_from_delta(chars_text: str, items_text: str) -> MorphologicalMatrix:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(chars_text)
        (p / "items").write_text(items_text)
        reader = DeltaReader(p)
        reader.read()
        return MatrixBuilder(reader.characters, reader.taxa).build()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def variable_matrix() -> MorphologicalMatrix:
    """3 taxa × 3 characters: one variable, one constant, one with missing."""
    df = pd.DataFrame(
        {
            "C1": [0.0, 1.0, 2.0],   # variable
            "C2": [1.0, 1.0, 1.0],   # constant
            "C3": [0.0, float("nan"), 1.0],  # variable + missing
        },
        index=["TaxA", "TaxB", "TaxC"],
    )
    return _make_matrix(df)


@pytest.fixture
def autamorphy_matrix() -> MorphologicalMatrix:
    """3 taxa × 3 characters.

    C1: 0,0,1  → shared derived state — NOT an autapomorphy (synapomorphy if 2-taxon)
    C2: 0,1,1  → state 0 unique to TaxA  (autapomorphy of TaxA for derived state)
    C3: 0,1,2  → each state unique (multi-autapomorphies)
    """
    df = pd.DataFrame(
        {
            "C1": [0.0, 0.0, 1.0],
            "C2": [0.0, 1.0, 1.0],
            "C3": [0.0, 1.0, 2.0],
        },
        index=["TaxA", "TaxB", "TaxC"],
    )
    return _make_matrix(df)


@pytest.fixture
def sparse_matrix() -> MorphologicalMatrix:
    """4 taxa × 3 characters with varying amounts of missing data."""
    df = pd.DataFrame(
        {
            "C1": [0.0, 1.0, float("nan"), float("nan")],   # 50% complete
            "C2": [0.0, float("nan"), float("nan"), float("nan")],  # 25% complete
            "C3": [0.0, 1.0, 0.0, 1.0],  # 100% complete
        },
        index=["T1", "T2", "T3", "T4"],
    )
    return _make_matrix(df)


# ---------------------------------------------------------------------------
# Tests: remove_constant_characters
# ---------------------------------------------------------------------------


class TestRemoveConstantCharacters:
    def test_removes_constant_column(self, variable_matrix):
        result = MatrixCleaner(variable_matrix).remove_constant_characters()
        assert "C2" not in result.df.columns

    def test_keeps_variable_column(self, variable_matrix):
        result = MatrixCleaner(variable_matrix).remove_constant_characters()
        assert "C1" in result.df.columns

    def test_keeps_column_with_missing_but_variable(self, variable_matrix):
        result = MatrixCleaner(variable_matrix).remove_constant_characters()
        assert "C3" in result.df.columns

    def test_output_shape_reduced(self, variable_matrix):
        result = MatrixCleaner(variable_matrix).remove_constant_characters()
        assert result.df.shape[1] < variable_matrix.df.shape[1]

    def test_all_constant_returns_empty(self):
        df = pd.DataFrame(
            {"C1": [1.0, 1.0], "C2": [2.0, 2.0]},
            index=["T1", "T2"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_constant_characters()
        assert result.df.shape[1] == 0

    def test_all_variable_nothing_removed(self):
        df = pd.DataFrame(
            {"C1": [0.0, 1.0], "C2": [1.0, 0.0]},
            index=["T1", "T2"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_constant_characters()
        assert result.df.shape[1] == 2

    def test_returns_new_matrix_not_same_object(self, variable_matrix):
        result = MatrixCleaner(variable_matrix).remove_constant_characters()
        assert result is not variable_matrix


# ---------------------------------------------------------------------------
# Tests: remove_autapomorphies
# ---------------------------------------------------------------------------


class TestRemoveAutapomorphies:
    def test_column_with_all_unique_states_removed(self, autamorphy_matrix):
        """C3 has three different states, each in exactly one taxon → autapomorphy."""
        result = MatrixCleaner(autamorphy_matrix).remove_autapomorphies()
        assert "C3" not in result.df.columns

    def test_synapomorphy_column_kept(self, autamorphy_matrix):
        """C2 has derived state 1 shared by TaxB and TaxC → synapomorphy, keep it."""
        result = MatrixCleaner(autamorphy_matrix).remove_autapomorphies()
        assert "C2" in result.df.columns

    def test_autapomorphy_single_derived_taxon_removed(self):
        """Character where exactly one taxon shows the derived state is autapomorphic."""
        df = pd.DataFrame(
            {"C1": [0.0, 0.0, 1.0]},
            index=["TaxA", "TaxB", "TaxC"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_autapomorphies()
        assert "C1" not in result.df.columns

    def test_synapomorphy_shared_by_two_taxa_kept(self):
        """Derived state shared by 2 taxa must not be removed."""
        df = pd.DataFrame(
            {"C1": [0.0, 1.0, 1.0]},
            index=["TaxA", "TaxB", "TaxC"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_autapomorphies()
        assert "C1" in result.df.columns

    def test_constant_ancestral_column_not_removed(self):
        """A constant column at state 0 has no derived states; must be preserved."""
        df = pd.DataFrame(
            {"C1": [0.0, 0.0, 0.0]},
            index=["T1", "T2", "T3"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_autapomorphies()
        assert "C1" in result.df.columns

    def test_does_not_remove_all_columns(self):
        """C1 (derived state shared by 2 taxa) must be kept; C2 is autapomorphic."""
        df = pd.DataFrame(
            {"C1": [0.0, 1.0, 1.0], "C2": [0.0, 1.0, 2.0]},
            index=["T1", "T2", "T3"],
        )
        result = MatrixCleaner(_make_matrix(df)).remove_autapomorphies()
        assert result.df.shape[1] >= 1
        assert "C1" in result.df.columns

    def test_returns_dataframe(self, autamorphy_matrix):
        result = MatrixCleaner(autamorphy_matrix).remove_autapomorphies()
        assert isinstance(result.df, pd.DataFrame)


# ---------------------------------------------------------------------------
# Tests: filter_by_completeness
# ---------------------------------------------------------------------------


class TestFilterByCompleteness:
    def test_removes_character_below_threshold(self, sparse_matrix):
        """C2 is only 25% complete; with min_char_completeness=0.5 it should go."""
        result = MatrixCleaner(sparse_matrix).filter_by_completeness(
            min_char_completeness=0.5, min_taxon_completeness=0.0
        )
        assert "C2" not in result.df.columns

    def test_keeps_fully_complete_character(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).filter_by_completeness(
            min_char_completeness=0.5, min_taxon_completeness=0.0
        )
        assert "C3" in result.df.columns

    def test_removes_taxon_below_threshold(self, sparse_matrix):
        """T3 and T4 have only C3 scored; at min_taxon_completeness=0.6 they go."""
        result = MatrixCleaner(sparse_matrix).filter_by_completeness(
            min_taxon_completeness=0.6, min_char_completeness=0.0
        )
        assert "T3" not in result.df.index
        assert "T4" not in result.df.index

    def test_threshold_zero_keeps_everything(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).filter_by_completeness(
            min_taxon_completeness=0.0, min_char_completeness=0.0
        )
        assert result.df.shape == sparse_matrix.df.shape

    def test_threshold_one_keeps_only_perfect(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).filter_by_completeness(
            min_taxon_completeness=1.0, min_char_completeness=1.0
        )
        # Only C3 is 100% complete; only taxa with C3 scored (all of them) remain.
        assert "C2" not in result.df.columns
        assert "C1" not in result.df.columns

    def test_returns_morphological_matrix(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).filter_by_completeness()
        assert isinstance(result, MorphologicalMatrix)


# ---------------------------------------------------------------------------
# Tests: impute_missing
# ---------------------------------------------------------------------------


class TestImputeMissing:
    def test_mode_imputation_fills_nan(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).impute_missing(strategy="mode")
        assert not result.df.isna().any().any()

    def test_zero_imputation_fills_with_zero(self, sparse_matrix):
        result = MatrixCleaner(sparse_matrix).impute_missing(strategy="zero")
        assert not result.df.isna().any().any()
        # Every cell that was NaN in the original should now be exactly 0.0
        was_nan = sparse_matrix.df.isna()
        for col in sparse_matrix.df.columns:
            for idx in sparse_matrix.df.index:
                if was_nan.loc[idx, col]:
                    assert result.df.loc[idx, col] == 0.0, (
                        f"Cell [{idx}, {col}] was NaN but is not 0 after zero imputation"
                    )

    def test_unknown_strategy_raises(self, sparse_matrix):
        with pytest.raises(ValueError, match="Unknown imputation strategy"):
            MatrixCleaner(sparse_matrix).impute_missing(strategy="median")

    def test_mode_uses_most_common_value(self):
        df = pd.DataFrame(
            {"C1": [1.0, 1.0, 0.0, float("nan")]},
            index=["T1", "T2", "T3", "T4"],
        )
        result = MatrixCleaner(_make_matrix(df)).impute_missing(strategy="mode")
        # Mode of C1 (ignoring NaN) is 1.0
        assert result.df.loc["T4", "C1"] == 1.0
