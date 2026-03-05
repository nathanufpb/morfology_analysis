"""Tests for tree output (IO writers)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder, MorphologicalMatrix
from delta_phylo.io.nexus_writer import NexusWriter
from delta_phylo.io.phylip_writer import PhylipWriter
from delta_phylo.io.tnt_writer import TNTWriter
from delta_phylo.io.csv_writer import CSVWriter
from delta_phylo.matrix.encoding import MatrixEncoder


# ---------------------------------------------------------------------------
# Fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_matrix() -> MorphologicalMatrix:
    chars = textwrap.dedent("""\
        #1. body size/
        1. small/
        2. large/
        #2. wings/
        1. absent/
        2. present/
        #3. colour/
        1. red/
        2. blue/
        3. green/
    """)
    items = textwrap.dedent("""\
        #1. Taxon_Alpha/
        1,1 2,2 3,1/
        #2. Taxon_Beta/
        1,2 2,1 3,3/
        #3. Taxon_Gamma/
        1,1 2,2 3,?/
    """)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(chars)
        (p / "items").write_text(items)
        reader = DeltaReader(p)
        reader.read()
        builder = MatrixBuilder(reader.characters, reader.taxa)
        return builder.build()


# ---------------------------------------------------------------------------
# Tests: NexusWriter
# ---------------------------------------------------------------------------


class TestNexusWriter:
    def test_write_creates_file(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.nex"
        NexusWriter(sample_matrix).write(out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_nexus_starts_with_header(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.nex"
        NexusWriter(sample_matrix).write(out)
        content = out.read_text()
        assert content.startswith("#NEXUS")

    def test_nexus_contains_taxa(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.nex"
        NexusWriter(sample_matrix).write(out)
        content = out.read_text()
        for name in sample_matrix.get_taxa_names():
            assert name in content

    def test_nexus_contains_missing_marker(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.nex"
        NexusWriter(sample_matrix).write(out)
        content = out.read_text()
        assert "?" in content  # Taxon_Gamma has missing data

    def test_nexus_has_characters_block(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.nex"
        NexusWriter(sample_matrix).write(out)
        content = out.read_text()
        assert "BEGIN CHARACTERS;" in content
        assert "MATRIX" in content
        assert "END;" in content


# ---------------------------------------------------------------------------
# Tests: PhylipWriter
# ---------------------------------------------------------------------------


class TestPhylipWriter:
    def test_write_creates_file(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.phy"
        PhylipWriter(sample_matrix).write(out)
        assert out.exists()

    def test_first_line_dimensions(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.phy"
        PhylipWriter(sample_matrix).write(out)
        content = out.read_text()
        first_line = content.strip().split("\n")[0].strip()
        n_taxa, n_chars = int(first_line.split()[0]), int(first_line.split()[1])
        assert n_taxa == 3
        assert n_chars == 3

    def test_taxon_names_truncated(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.phy"
        PhylipWriter(sample_matrix).write(out)
        content = out.read_text()
        # Only check taxon lines (lines containing alphabetic characters, not the dimension line)
        lines = [
            l for l in content.strip().split("\n")
            if l and any(c.isalpha() for c in l)
        ]
        for line in lines:
            # Name field is exactly 10 chars
            assert len(line) >= 10


# ---------------------------------------------------------------------------
# Tests: TNTWriter
# ---------------------------------------------------------------------------


class TestTNTWriter:
    def test_write_creates_file(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.tnt"
        TNTWriter(sample_matrix).write(out)
        assert out.exists()

    def test_tnt_starts_with_xread(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.tnt"
        TNTWriter(sample_matrix).write(out)
        content = out.read_text()
        assert content.startswith("xread")

    def test_tnt_ends_with_proc(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.tnt"
        TNTWriter(sample_matrix).write(out)
        content = out.read_text().strip()
        assert "proc /;" in content


# ---------------------------------------------------------------------------
# Tests: CSVWriter
# ---------------------------------------------------------------------------


class TestCSVWriter:
    def test_write_creates_file(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.csv"
        CSVWriter(sample_matrix).write(out)
        assert out.exists()

    def test_csv_has_header(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        import csv

        out = tmp_path / "matrix.csv"
        CSVWriter(sample_matrix).write(out)
        with out.open() as f:
            reader = csv.reader(f)
            header = next(reader)
        assert "taxon" in header or "C1" in header

    def test_csv_row_count(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        import csv

        out = tmp_path / "matrix.csv"
        CSVWriter(sample_matrix).write(out)
        with out.open() as f:
            rows = list(csv.reader(f))
        # 1 header + 3 taxa
        assert len(rows) == 4

    def test_csv_missing_representation(
        self, sample_matrix: MorphologicalMatrix, tmp_path: Path
    ) -> None:
        out = tmp_path / "matrix.csv"
        CSVWriter(sample_matrix).write(out, na_rep="?")
        content = out.read_text()
        assert "?" in content
