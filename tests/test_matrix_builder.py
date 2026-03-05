"""Tests for the matrix builder and cleaning modules."""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder, MorphologicalMatrix
from delta_phylo.matrix.matrix_cleaning import MatrixCleaner
from delta_phylo.matrix.encoding import MatrixEncoder


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_matrix() -> MorphologicalMatrix:
    """Build a simple 3-taxon × 4-character matrix for testing."""
    delta_text_chars = textwrap.dedent("""\
        #1. body size/
        1. small/
        2. medium/
        3. large/
        #2. wings/
        1. absent/
        2. present/
        #3. legs/
        1. six/
        2. eight/
        3. many/
        #4. colour/
        1. red/
        2. blue/
    """)
    delta_text_items = textwrap.dedent("""\
        #1. TaxonA/
        1,1 2,2 3,1 4,1/
        #2. TaxonB/
        1,2 2,1 3,2 4,2/
        #3. TaxonC/
        1,3 2,2 3,3 4,?/
    """)
    import tempfile
    import os

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(delta_text_chars)
        (p / "items").write_text(delta_text_items)
        reader = DeltaReader(p)
        reader.read()
        builder = MatrixBuilder(reader.characters, reader.taxa)
        return builder.build()


@pytest.fixture
def constant_matrix() -> MorphologicalMatrix:
    """Matrix where one character is constant across all taxa."""
    delta_text_chars = textwrap.dedent("""\
        #1. variable/
        1. state_a/
        2. state_b/
        #2. constant/
        1. only_state/
        2. placeholder/
    """)
    delta_text_items = textwrap.dedent("""\
        #1. TaxonA/
        1,1 2,1/
        #2. TaxonB/
        1,2 2,1/
        #3. TaxonC/
        1,1 2,1/
    """)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(delta_text_chars)
        (p / "items").write_text(delta_text_items)
        reader = DeltaReader(p)
        reader.read()
        builder = MatrixBuilder(reader.characters, reader.taxa)
        return builder.build()


# ---------------------------------------------------------------------------
# Tests: MatrixBuilder
# ---------------------------------------------------------------------------


class TestMatrixBuilder:
    def test_shape(self, simple_matrix: MorphologicalMatrix) -> None:
        assert simple_matrix.shape == (3, 4)

    def test_taxa_names_in_index(self, simple_matrix: MorphologicalMatrix) -> None:
        names = simple_matrix.get_taxa_names()
        assert "TaxonA" in names
        assert "TaxonB" in names
        assert "TaxonC" in names

    def test_column_names_start_with_c(self, simple_matrix: MorphologicalMatrix) -> None:
        cols = simple_matrix.get_character_names()
        for col in cols:
            assert col.startswith("C"), f"Expected column to start with 'C', got {col}"

    def test_missing_data_is_nan(self, simple_matrix: MorphologicalMatrix) -> None:
        # TaxonC character 4 is missing (?)
        val = simple_matrix.df.loc["TaxonC", "C4"]
        assert pd.isna(val)

    def test_state_conversion(self, simple_matrix: MorphologicalMatrix) -> None:
        # DELTA state 1 → 0-based index 0
        val = simple_matrix.df.loc["TaxonA", "C1"]
        assert val == 0.0

    def test_to_numpy_shape(self, simple_matrix: MorphologicalMatrix) -> None:
        arr = simple_matrix.to_numpy()
        assert arr.shape == (3, 4)

    def test_to_numpy_missing_sentinel(self, simple_matrix: MorphologicalMatrix) -> None:
        arr = simple_matrix.to_numpy(missing_value=-1)
        # TaxonC's character 4 should be -1
        taxon_names = simple_matrix.get_taxa_names()
        taxon_c_idx = taxon_names.index("TaxonC")
        assert arr[taxon_c_idx, 3] == -1

    def test_polymorphic_first_strategy(self) -> None:
        """Test that polymorphic states resolve to the lowest state with 'first' strategy."""
        import tempfile

        delta_chars = "#1. colour/\n1. red/\n2. blue/\n3. green/\n"
        delta_items = "#1. TaxonPoly/\n1,2/3/\n"
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "characters").write_text(delta_chars)
            (p / "items").write_text(delta_items)
            reader = DeltaReader(p)
            reader.read()
            builder = MatrixBuilder(
                reader.characters, reader.taxa, polymorphic_strategy="first"
            )
            matrix = builder.build()
        # States 2 and 3 in DELTA → 0-based indices 1 and 2 → min = 1
        val = matrix.df.iloc[0, 0]
        assert val == 1.0

    def test_polymorphic_missing_strategy(self) -> None:
        """Test that polymorphic states become NaN with 'missing' strategy."""
        import tempfile

        delta_chars = "#1. colour/\n1. red/\n2. blue/\n3. green/\n"
        delta_items = "#1. TaxonPoly/\n1,2/3/\n"
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "characters").write_text(delta_chars)
            (p / "items").write_text(delta_items)
            reader = DeltaReader(p)
            reader.read()
            builder = MatrixBuilder(
                reader.characters, reader.taxa, polymorphic_strategy="missing"
            )
            matrix = builder.build()
        val = matrix.df.iloc[0, 0]
        assert pd.isna(val)


# ---------------------------------------------------------------------------
# Tests: MatrixCleaner
# ---------------------------------------------------------------------------


class TestMatrixCleaner:
    def test_remove_constant_characters(
        self, constant_matrix: MorphologicalMatrix
    ) -> None:
        cleaner = MatrixCleaner(constant_matrix)
        cleaned = cleaner.remove_constant_characters()
        # C2 is constant (all state 0), should be removed
        assert cleaned.shape[1] == 1
        assert "C1" in cleaned.get_character_names()
        assert "C2" not in cleaned.get_character_names()

    def test_filter_by_completeness_all_kept(
        self, simple_matrix: MorphologicalMatrix
    ) -> None:
        cleaner = MatrixCleaner(simple_matrix)
        # With low thresholds everything should be kept
        filtered = cleaner.filter_by_completeness(
            min_taxon_completeness=0.0, min_char_completeness=0.0
        )
        assert filtered.shape == simple_matrix.shape

    def test_impute_missing_mode(self, simple_matrix: MorphologicalMatrix) -> None:
        cleaner = MatrixCleaner(simple_matrix)
        imputed = cleaner.impute_missing(strategy="mode")
        # After imputation, no NaN should remain
        assert not imputed.df.isna().any().any()

    def test_impute_missing_zero(self, simple_matrix: MorphologicalMatrix) -> None:
        cleaner = MatrixCleaner(simple_matrix)
        imputed = cleaner.impute_missing(strategy="zero")
        assert not imputed.df.isna().any().any()
        # Previously NaN cell should now be 0
        assert imputed.df.loc["TaxonC", "C4"] == 0.0


# ---------------------------------------------------------------------------
# Tests: MatrixEncoder
# ---------------------------------------------------------------------------


class TestMatrixEncoder:
    def test_symbol_string_no_missing(self, simple_matrix: MorphologicalMatrix) -> None:
        encoder = MatrixEncoder(simple_matrix)
        sym = encoder.symbol_string("TaxonA")
        assert isinstance(sym, str)
        assert "?" not in sym

    def test_symbol_string_with_missing(
        self, simple_matrix: MorphologicalMatrix
    ) -> None:
        encoder = MatrixEncoder(simple_matrix)
        sym = encoder.symbol_string("TaxonC")
        assert "?" in sym

    def test_symbol_string_unknown_taxon(
        self, simple_matrix: MorphologicalMatrix
    ) -> None:
        encoder = MatrixEncoder(simple_matrix)
        with pytest.raises(KeyError):
            encoder.symbol_string("DoesNotExist")

    def test_to_binary_splits_multistate(
        self, simple_matrix: MorphologicalMatrix
    ) -> None:
        encoder = MatrixEncoder(simple_matrix)
        binary = encoder.to_binary()
        # Original C1 has 3 states → should become 3 binary columns
        n_original = simple_matrix.shape[1]
        assert binary.shape[1] >= n_original
