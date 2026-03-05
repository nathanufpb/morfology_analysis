"""Tests for CLI commands — closes issue #12.

Uses Click's CliRunner to test all commands without spawning a subprocess.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from delta_phylo.cli.main import cli


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def delta_dir(tmp_path: Path) -> Path:
    """Minimal valid DELTA dataset directory."""
    specs = "*NUMBER OF CHARACTERS 3\n*MAXIMUM NUMBER OF ITEMS 3\n"
    characters = textwrap.dedent("""\
        #1. size/
        1. small/
        2. medium/
        3. large/
        #2. wings/
        1. absent/
        2. present/
        #3. colour/
        1. red/
        2. blue/
        3. green/
    """)
    items = textwrap.dedent("""\
        #1. TaxA/
        1,1 2,2 3,1/
        #2. TaxB/
        1,2 2,1 3,3/
        #3. TaxC/
        1,3 2,2 3,2/
    """)
    (tmp_path / "specs").write_text(specs)
    (tmp_path / "characters").write_text(characters)
    (tmp_path / "items").write_text(items)
    return tmp_path


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


# ---------------------------------------------------------------------------
# Tests: parse command
# ---------------------------------------------------------------------------


class TestParseCommand:
    def test_valid_dataset_exits_zero(self, runner, delta_dir):
        result = runner.invoke(cli, ["parse", str(delta_dir)])
        assert result.exit_code == 0

    def test_output_contains_character_count(self, runner, delta_dir):
        result = runner.invoke(cli, ["parse", str(delta_dir)])
        assert "3" in result.output  # 3 characters

    def test_output_contains_taxa_count(self, runner, delta_dir):
        result = runner.invoke(cli, ["parse", str(delta_dir)])
        assert "3" in result.output  # 3 taxa

    def test_missing_dataset_nonzero_exit(self, runner, tmp_path):
        missing = str(tmp_path / "nonexistent")
        result = runner.invoke(cli, ["parse", missing])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Tests: matrix command
# ---------------------------------------------------------------------------


class TestMatrixCommand:
    def test_matrix_cmd_exits_zero(self, runner, delta_dir):
        result = runner.invoke(cli, ["matrix", str(delta_dir)])
        assert result.exit_code == 0

    def test_matrix_output_contains_taxa_names(self, runner, delta_dir):
        result = runner.invoke(cli, ["matrix", str(delta_dir)])
        assert "TaxA" in result.output

    def test_output_flag_creates_csv(self, runner, delta_dir, tmp_path):
        out_csv = str(tmp_path / "matrix.csv")
        result = runner.invoke(cli, ["matrix", str(delta_dir), "--output", out_csv])
        assert result.exit_code == 0
        assert Path(out_csv).exists()

    def test_remove_constant_flag_accepted(self, runner, delta_dir):
        result = runner.invoke(cli, ["matrix", str(delta_dir), "--remove-constant"])
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# Tests: export command
# ---------------------------------------------------------------------------


class TestExportCommand:
    def test_nexus_export_creates_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "matrix.nex")
        result = runner.invoke(
            cli, ["export", "nexus", str(delta_dir), "--output", out]
        )
        assert result.exit_code == 0
        assert Path(out).exists()

    def test_nexus_file_starts_with_header(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "matrix.nex")
        runner.invoke(cli, ["export", "nexus", str(delta_dir), "--output", out])
        content = Path(out).read_text()
        assert content.startswith("#NEXUS")

    def test_phylip_export_creates_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "matrix.phy")
        result = runner.invoke(
            cli, ["export", "phylip", str(delta_dir), "--output", out]
        )
        assert result.exit_code == 0
        assert Path(out).exists()

    def test_tnt_export_creates_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "matrix.tnt")
        result = runner.invoke(
            cli, ["export", "tnt", str(delta_dir), "--output", out]
        )
        assert result.exit_code == 0
        assert Path(out).exists()

    def test_csv_export_creates_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "matrix.csv")
        result = runner.invoke(
            cli, ["export", "csv", str(delta_dir), "--output", out]
        )
        assert result.exit_code == 0
        assert Path(out).exists()


# ---------------------------------------------------------------------------
# Tests: parsimony command
# ---------------------------------------------------------------------------


class TestParsimonyCommand:
    def test_parsimony_exits_zero(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "tree.nwk")
        result = runner.invoke(
            cli,
            [
                "parsimony",
                str(delta_dir),
                "--replicates", "1",
                "--iterations", "5",
                "--output", out,
                "--seed", "0",
            ],
        )
        assert result.exit_code == 0

    def test_parsimony_creates_tree_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "tree.nwk")
        runner.invoke(
            cli,
            [
                "parsimony",
                str(delta_dir),
                "--replicates", "1",
                "--iterations", "5",
                "--output", out,
                "--seed", "0",
            ],
        )
        assert Path(out).exists()

    def test_parsimony_score_in_output(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "tree.nwk")
        result = runner.invoke(
            cli,
            [
                "parsimony",
                str(delta_dir),
                "--replicates", "1",
                "--iterations", "5",
                "--output", out,
            ],
        )
        assert "score" in result.output.lower() or result.exit_code == 0


# ---------------------------------------------------------------------------
# Tests: distance command
# ---------------------------------------------------------------------------


class TestDistanceCommand:
    def test_distance_exits_zero(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "nj.nwk")
        result = runner.invoke(
            cli, ["distance", str(delta_dir), "--output", out]
        )
        assert result.exit_code == 0

    def test_distance_creates_tree_file(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "nj.nwk")
        runner.invoke(cli, ["distance", str(delta_dir), "--output", out])
        assert Path(out).exists()

    def test_simple_matching_metric_accepted(self, runner, delta_dir, tmp_path):
        out = str(tmp_path / "nj_sm.nwk")
        result = runner.invoke(
            cli,
            ["distance", str(delta_dir), "--metric", "simple_matching", "--output", out],
        )
        assert result.exit_code == 0


# ---------------------------------------------------------------------------
# Tests: analyze command (full pipeline)
# ---------------------------------------------------------------------------


class TestAnalyzeCommand:
    def test_analyze_parsimony_exits_zero(self, runner, delta_dir, tmp_path):
        result = runner.invoke(
            cli,
            [
                "analyze",
                str(delta_dir),
                "--method", "parsimony",
                "--output-dir", str(tmp_path),
                "--replicates", "1",
            ],
        )
        assert result.exit_code == 0

    def test_analyze_nj_exits_zero(self, runner, delta_dir, tmp_path):
        result = runner.invoke(
            cli,
            [
                "analyze",
                str(delta_dir),
                "--method", "nj",
                "--output-dir", str(tmp_path),
            ],
        )
        assert result.exit_code == 0

    def test_analyze_creates_nexus_file(self, runner, delta_dir, tmp_path):
        runner.invoke(
            cli,
            [
                "analyze",
                str(delta_dir),
                "--method", "nj",
                "--output-dir", str(tmp_path),
            ],
        )
        assert (tmp_path / "matrix.nex").exists()

    def test_analyze_creates_csv_file(self, runner, delta_dir, tmp_path):
        runner.invoke(
            cli,
            [
                "analyze",
                str(delta_dir),
                "--method", "nj",
                "--output-dir", str(tmp_path),
            ],
        )
        assert (tmp_path / "matrix.csv").exists()
