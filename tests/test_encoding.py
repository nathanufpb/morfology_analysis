"""Tests for MatrixEncoder — closes issue #8.

Covers:
  - to_binary(): multistate → one-hot expansion
  - symbol_string(): correct symbols, missing data, large states
  - _encode_state(): static helper
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix, MatrixBuilder
from delta_phylo.matrix.encoding import MatrixEncoder
from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.parser.characters import Character, CharacterType


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _matrix_from_delta(chars_text: str, items_text: str) -> MorphologicalMatrix:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(chars_text)
        (p / "items").write_text(items_text)
        reader = DeltaReader(p)
        reader.read()
        return MatrixBuilder(reader.characters, reader.taxa).build()


def _make_matrix(df: pd.DataFrame, characters=None) -> MorphologicalMatrix:
    return MorphologicalMatrix(df=df, characters=characters or [], taxa=[])


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def binary_matrix() -> MorphologicalMatrix:
    """2 taxa × 2 binary characters."""
    chars = textwrap.dedent("""\
        #1. wings/ 1. absent/ 2. present/
        #2. colour/ 1. red/ 2. blue/
    """)
    items = textwrap.dedent("""\
        #1. TaxA/ 1,1 2,2/
        #2. TaxB/ 1,2 2,1/
    """)
    return _matrix_from_delta(chars, items)


@pytest.fixture
def multistate_matrix() -> MorphologicalMatrix:
    """2 taxa × 1 three-state character."""
    chars = textwrap.dedent("""\
        #1. size/ 1. small/ 2. medium/ 3. large/
    """)
    items = textwrap.dedent("""\
        #1. TaxA/ 1,1/
        #2. TaxB/ 1,3/
    """)
    return _matrix_from_delta(chars, items)


@pytest.fixture
def missing_matrix() -> MorphologicalMatrix:
    """2 taxa × 1 character where one taxon has missing data."""
    chars = textwrap.dedent("""\
        #1. wings/ 1. absent/ 2. present/
    """)
    items = textwrap.dedent("""\
        #1. TaxA/ 1,1/
        #2. TaxB/ 1,?/
    """)
    return _matrix_from_delta(chars, items)


# ---------------------------------------------------------------------------
# Tests: symbol_string
# ---------------------------------------------------------------------------


class TestSymbolString:
    def test_correct_state_symbol_for_known_taxon(self, binary_matrix):
        encoder = MatrixEncoder(binary_matrix)
        sym = encoder.symbol_string("TaxA")
        # TaxA: wings=1 (state index 0), colour=2 (state index 1)
        assert sym[0] == "0"   # state 0
        assert sym[1] == "1"   # state 1

    def test_missing_data_encoded_as_question_mark(self, missing_matrix):
        encoder = MatrixEncoder(missing_matrix)
        sym = encoder.symbol_string("TaxB")
        assert sym == "?"

    def test_unknown_taxon_raises_key_error(self, binary_matrix):
        encoder = MatrixEncoder(binary_matrix)
        with pytest.raises(KeyError):
            encoder.symbol_string("NonExistentTaxon")

    def test_symbol_string_length_equals_num_characters(self, binary_matrix):
        encoder = MatrixEncoder(binary_matrix)
        sym = encoder.symbol_string("TaxA")
        assert len(sym) == len(binary_matrix.get_character_names())


class TestEncodeState:
    def test_single_digit_states(self):
        for i in range(10):
            assert MatrixEncoder._encode_state(i) == str(i)

    def test_state_ten_is_capital_a(self):
        assert MatrixEncoder._encode_state(10) == "A"

    def test_state_eleven_is_capital_b(self):
        assert MatrixEncoder._encode_state(11) == "B"

    def test_state_35_is_capital_z(self):
        assert MatrixEncoder._encode_state(35) == "Z"

    def test_negative_state_encodes_as_question_mark(self):
        assert MatrixEncoder._encode_state(-1) == "?"


# ---------------------------------------------------------------------------
# Tests: to_binary
# ---------------------------------------------------------------------------


class TestToBinary:
    def test_binary_columns_unchanged(self, binary_matrix):
        """Binary characters should not be expanded."""
        encoder = MatrixEncoder(binary_matrix)
        result = encoder.to_binary()
        # Original columns should still be present (not expanded)
        assert result.df.shape[0] == binary_matrix.df.shape[0]
        # Total columns = same (all are binary)
        assert result.df.shape[1] == binary_matrix.df.shape[1]

    def test_three_state_character_produces_three_columns(self, multistate_matrix):
        """A 3-state character should expand into 3 binary sub-columns."""
        encoder = MatrixEncoder(multistate_matrix)
        result = encoder.to_binary()
        # Original had 1 column (C1 with max_state=2); should expand to 3
        assert result.df.shape[1] == 3

    def test_binary_values_only_zero_or_one(self, multistate_matrix):
        encoder = MatrixEncoder(multistate_matrix)
        result = encoder.to_binary()
        valid_values = result.df.dropna().values.flatten()
        assert set(valid_values).issubset({0.0, 1.0})

    def test_nan_propagated_in_binary_expansion(self, missing_matrix):
        encoder = MatrixEncoder(missing_matrix)
        result = encoder.to_binary()
        # TaxB had missing data; all its binary columns should be NaN
        assert result.df.loc["TaxB"].isna().all()

    def test_returns_morphological_matrix(self, multistate_matrix):
        result = MatrixEncoder(multistate_matrix).to_binary()
        assert isinstance(result, MorphologicalMatrix)

    def test_row_count_preserved(self, multistate_matrix):
        result = MatrixEncoder(multistate_matrix).to_binary()
        assert result.df.shape[0] == multistate_matrix.df.shape[0]
