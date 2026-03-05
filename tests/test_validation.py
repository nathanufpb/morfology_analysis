"""Tests for validation utilities — closes issue #10.

Covers:
  - validate_delta_directory()
  - validate_matrix()
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from delta_phylo.utils.validation import validate_delta_directory, validate_matrix
from delta_phylo.matrix.matrix_builder import MorphologicalMatrix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_matrix(df: pd.DataFrame) -> MorphologicalMatrix:
    return MorphologicalMatrix(df=df, characters=[], taxa=[])


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def valid_delta_dir(tmp_path: Path) -> Path:
    """A minimal valid DELTA directory."""
    (tmp_path / "items").write_text("#1. TaxA/ 1,1/\n")
    return tmp_path


# ---------------------------------------------------------------------------
# Tests: validate_delta_directory
# ---------------------------------------------------------------------------


class TestValidateDeltaDirectory:
    def test_valid_directory_returns_path(self, valid_delta_dir):
        result = validate_delta_directory(valid_delta_dir)
        assert isinstance(result, Path)
        assert result.exists()

    def test_returned_path_is_resolved(self, valid_delta_dir):
        relative = valid_delta_dir  # already absolute in tmp_path
        result = validate_delta_directory(relative)
        assert result == result.resolve()

    def test_nonexistent_path_raises_file_not_found(self, tmp_path):
        missing = tmp_path / "does_not_exist"
        with pytest.raises(FileNotFoundError):
            validate_delta_directory(missing)

    def test_file_path_raises_not_a_directory(self, tmp_path):
        a_file = tmp_path / "somefile.txt"
        a_file.write_text("hello")
        with pytest.raises(NotADirectoryError):
            validate_delta_directory(a_file)

    def test_empty_directory_raises_value_error(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(ValueError, match="No recognised DELTA files"):
            validate_delta_directory(empty_dir)

    def test_dir_without_recognised_files_raises_value_error(self, tmp_path):
        d = tmp_path / "norecognised"
        d.mkdir()
        (d / "random.txt").write_text("data")
        with pytest.raises(ValueError):
            validate_delta_directory(d)

    def test_accepts_string_path(self, valid_delta_dir):
        result = validate_delta_directory(str(valid_delta_dir))
        assert isinstance(result, Path)

    def test_characters_file_alone_is_valid(self, tmp_path):
        (tmp_path / "characters").write_text("#1. size/ 1. small/\n")
        # Should not raise — 'characters' is a recognised filename
        result = validate_delta_directory(tmp_path)
        assert isinstance(result, Path)


# ---------------------------------------------------------------------------
# Tests: validate_matrix
# ---------------------------------------------------------------------------


class TestValidateMatrix:
    def test_valid_matrix_does_not_raise(self):
        df = pd.DataFrame(
            {"C1": [0.0, 1.0], "C2": [1.0, 0.0]},
            index=["T1", "T2"],
        )
        # Should not raise
        validate_matrix(_make_matrix(df))

    def test_single_taxon_raises_value_error(self):
        df = pd.DataFrame({"C1": [0.0]}, index=["T1"])
        with pytest.raises(ValueError, match="taxa"):
            validate_matrix(_make_matrix(df), min_taxa=2)

    def test_zero_characters_raises_value_error(self):
        df = pd.DataFrame(index=["T1", "T2"])
        with pytest.raises(ValueError, match="characters"):
            validate_matrix(_make_matrix(df))

    def test_min_taxa_respects_parameter(self):
        df = pd.DataFrame(
            {"C1": [0.0, 1.0, 2.0]},
            index=["T1", "T2", "T3"],
        )
        # min_taxa=3 should pass; min_taxa=4 should fail
        validate_matrix(_make_matrix(df), min_taxa=3)
        with pytest.raises(ValueError):
            validate_matrix(_make_matrix(df), min_taxa=4)

    def test_high_missing_data_logs_warning(self, caplog):
        """A matrix with >90% missing data should emit a WARNING log message."""
        # 10 taxa × 10 characters; only 1 cell is scored → 99% missing > 90%
        data = {f"C{i}": [float("nan")] * 10 for i in range(1, 11)}
        data["C1"][0] = 1.0  # single non-NaN value → 99/100 = 99% missing
        df = pd.DataFrame(data, index=[f"T{i}" for i in range(10)])
        with caplog.at_level(logging.WARNING, logger="delta_phylo.utils.validation"):
            validate_matrix(_make_matrix(df))
        assert any(
            "missing" in msg.lower() or "90" in msg or "%" in msg
            for msg in caplog.messages
        )

    def test_fully_missing_matrix_raises(self):
        """A completely empty matrix (0 chars) raises ValueError."""
        df = pd.DataFrame(index=["T1", "T2"])
        with pytest.raises(ValueError):
            validate_matrix(_make_matrix(df))
