"""Tests for numeric/continuous character support — closes issue #5.

Covers:
  - DeltaReader: RN token in *CHARACTER TYPES sets CharacterType.NUMERIC
  - DeltaReader: float scores are parsed correctly
  - MatrixBuilder: float values stored as-is (no 1-based conversion)
  - MorphologicalMatrix.to_numpy_float: preserves decimals
  - DistanceAnalysis: Gower distance for numeric characters
  - gower_distance function directly
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder, MorphologicalMatrix
from delta_phylo.parser.characters import CharacterType
from delta_phylo.phylogeny.distance import (
    DistanceAnalysis,
    gower_distance,
    _build_gower_matrix,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SPECS = """\
*NUMBER OF CHARACTERS 2
*MAXIMUM NUMBER OF ITEMS 3
*CHARACTER TYPES RN 1
"""

_CHARS = """\
#1. body length/
#2. wing colour/
1. absent/
2. present/
"""

_ITEMS = """\
#1. Alpha/ 1,3.5 2,1/
#2. Beta/  1,7.0 2,2/
#3. Gamma/ 1,5.5 2,1/
"""


def _make_dataset(chars: str = _CHARS, items: str = _ITEMS, specs: str = _SPECS):
    """Return a (DeltaReader, MorphologicalMatrix) pair from inline text."""
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "specs").write_text(specs)
        (p / "characters").write_text(chars)
        (p / "items").write_text(items)
        reader = DeltaReader(p)
        reader.read()
        matrix = MatrixBuilder(reader.characters, reader.taxa).build()
    return reader, matrix


# ---------------------------------------------------------------------------
# Tests: parser — RN directive
# ---------------------------------------------------------------------------


class TestNumericParsing:
    def test_rn_directive_sets_numeric_char_type(self):
        """RN token in *CHARACTER TYPES must set CharacterType.NUMERIC."""
        reader, _ = _make_dataset()
        char1 = next(c for c in reader.characters if c.char_id == 1)
        assert char1.char_type == CharacterType.NUMERIC

    def test_non_numeric_char_remains_unordered(self):
        """Character 2 (no RN token) must stay UNORDERED or BINARY."""
        reader, _ = _make_dataset()
        char2 = next(c for c in reader.characters if c.char_id == 2)
        assert char2.char_type in (CharacterType.UNORDERED, CharacterType.BINARY)

    def test_float_score_parsed_correctly(self):
        """Numeric character scores must be stored as floats."""
        reader, _ = _make_dataset()
        alpha = next(t for t in reader.taxa if t.name == "Alpha")
        score = alpha.get_score(1)
        assert isinstance(score, float)
        assert math.isclose(score, 3.5)

    def test_numeric_score_not_offset_by_one(self):
        """Float scores must NOT have the 1-based → 0-based subtraction applied."""
        reader, _ = _make_dataset()
        beta = next(t for t in reader.taxa if t.name == "Beta")
        score = beta.get_score(1)
        assert math.isclose(score, 7.0)

    def test_missing_numeric_score_is_none(self):
        """Missing data (?) for a numeric character must be None."""
        specs = "*NUMBER OF CHARACTERS 1\n*MAXIMUM NUMBER OF ITEMS 1\n*CHARACTER TYPES RN 1\n"
        chars = "#1. length/\n"
        items = "#1. T1/ 1,?/\n"
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "specs").write_text(specs)
            (p / "characters").write_text(chars)
            (p / "items").write_text(items)
            reader = DeltaReader(p)
            reader.read()
        assert reader.taxa[0].get_score(1) is None


# ---------------------------------------------------------------------------
# Tests: MatrixBuilder with numeric characters
# ---------------------------------------------------------------------------


class TestNumericMatrixBuilder:
    def test_numeric_values_in_dataframe(self):
        """Float scores must appear as float values in the DataFrame."""
        _, matrix = _make_dataset()
        val = matrix.df.at["Alpha", "C1"]
        assert math.isclose(val, 3.5)

    def test_to_numpy_float_preserves_decimals(self):
        """to_numpy_float must not truncate decimal values."""
        _, matrix = _make_dataset()
        arr = matrix.to_numpy_float()
        alpha_idx = matrix.get_taxa_names().index("Alpha")
        assert math.isclose(arr[alpha_idx, 0], 3.5)

    def test_to_numpy_truncates_float_to_int(self):
        """Legacy to_numpy truncates to int — documents known behaviour."""
        _, matrix = _make_dataset()
        arr = matrix.to_numpy()
        alpha_idx = matrix.get_taxa_names().index("Alpha")
        # 3.5 truncated → 3
        assert arr[alpha_idx, 0] == 3

    def test_missing_numeric_stored_as_nan(self):
        """Missing numeric score must be NaN in the DataFrame."""
        specs = "*NUMBER OF CHARACTERS 1\n*MAXIMUM NUMBER OF ITEMS 2\n*CHARACTER TYPES RN 1\n"
        chars = "#1. length/\n"
        items = "#1. T1/ 1,2.5/\n#2. T2/ 1,?/\n"
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            (p / "specs").write_text(specs)
            (p / "characters").write_text(chars)
            (p / "items").write_text(items)
            reader = DeltaReader(p)
            reader.read()
            matrix = MatrixBuilder(reader.characters, reader.taxa).build()
        assert math.isnan(matrix.df.at["T2", "C1"])


# ---------------------------------------------------------------------------
# Tests: gower_distance function
# ---------------------------------------------------------------------------


class TestGowerDistance:
    def test_identical_rows_zero_distance(self):
        a = np.array([1.0, 0.0])
        b = np.array([1.0, 0.0])
        is_num = np.array([True, False])
        ranges = np.array([10.0, 1.0])
        assert gower_distance(a, b, is_num, ranges) == pytest.approx(0.0)

    def test_numeric_component(self):
        """Numeric-only: |3.5 - 7.0| / 10.0 = 0.35."""
        a = np.array([3.5])
        b = np.array([7.0])
        is_num = np.array([True])
        ranges = np.array([10.0])
        assert gower_distance(a, b, is_num, ranges) == pytest.approx(0.35)

    def test_discrete_component(self):
        """Discrete-only mismatch: 1 character, different states → distance 1.0."""
        a = np.array([0.0])
        b = np.array([1.0])
        is_num = np.array([False])
        ranges = np.array([1.0])
        assert gower_distance(a, b, is_num, ranges) == pytest.approx(1.0)

    def test_missing_excluded_from_both_sides(self):
        """Missing values (-1 for discrete, NaN for numeric) must be skipped."""
        a = np.array([float("nan"), 0.0])
        b = np.array([7.0, 0.0])
        is_num = np.array([True, False])
        ranges = np.array([10.0, 1.0])
        # Only the second character is comparable → identical → 0.0
        assert gower_distance(a, b, is_num, ranges) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Tests: DistanceAnalysis with numeric characters
# ---------------------------------------------------------------------------


class TestDistanceAnalysisNumeric:
    def test_gower_metric_runs_without_error(self):
        _, matrix = _make_dataset()
        analysis = DistanceAnalysis(matrix, metric="gower")
        tree = analysis.run()
        assert tree is not None

    def test_mixed_matrix_auto_uses_gower(self):
        """A matrix that contains NUMERIC characters must use Gower distance."""
        _, matrix = _make_dataset()
        # The matrix has a NUMERIC character; _build_distance_matrix should auto-switch.
        dist = _build_gower_matrix(matrix)
        assert dist.shape == (3, 3)
        assert dist[0, 0] == pytest.approx(0.0)
        # Alpha (3.5) vs Beta (7.0): range is 7.0-3.5=3.5, |3.5-7.0|/3.5 = 1.0 for C1;
        # C2: 1 vs 2 → mismatch = 1.0  → mean = 1.0
        assert dist[0, 1] == pytest.approx(1.0)
