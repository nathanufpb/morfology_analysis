"""Integration tests using the bundled example DELTA dataset — closes issue #13.

Exercises the full pipeline:
  parse → matrix build → cleaning → distance NJ → parsimony → NEXUS export.

All tests use the real files in examples/example_delta_dataset/.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Phylo.BaseTree import Tree

# Path to the bundled example dataset (relative to this file's location)
EXAMPLE_DIR = Path(__file__).parent.parent / "examples" / "example_delta_dataset"


# ---------------------------------------------------------------------------
# Skip guard: skip all tests if the example dataset is missing
# ---------------------------------------------------------------------------


pytestmark = pytest.mark.skipif(
    not EXAMPLE_DIR.exists(),
    reason=f"Example dataset not found at {EXAMPLE_DIR}",
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def reader():
    from delta_phylo.parser.delta_reader import DeltaReader

    r = DeltaReader(EXAMPLE_DIR)
    r.read()
    return r


@pytest.fixture(scope="module")
def matrix(reader):
    from delta_phylo.matrix.matrix_builder import MatrixBuilder

    return MatrixBuilder(reader.characters, reader.taxa).build()


# ---------------------------------------------------------------------------
# Tests: parsing
# ---------------------------------------------------------------------------


class TestIntegrationParsing:
    def test_parses_correct_number_of_characters(self, reader):
        # specs declares 20 characters
        assert len(reader.characters) == 20

    def test_parses_correct_number_of_taxa(self, reader):
        # specs declares 10 items
        assert len(reader.taxa) == 10

    def test_all_taxa_have_names(self, reader):
        for taxon in reader.taxa:
            assert taxon.name, f"Taxon {taxon.taxon_id} has no name"

    def test_all_characters_have_names(self, reader):
        for char in reader.characters:
            assert char.name, f"Character {char.char_id} has no name"

    def test_taxon_ids_sequential(self, reader):
        ids = [t.taxon_id for t in reader.taxa]
        assert ids == list(range(1, len(ids) + 1))

    def test_character_ids_sequential(self, reader):
        ids = [c.char_id for c in reader.characters]
        assert ids == list(range(1, len(ids) + 1))


# ---------------------------------------------------------------------------
# Tests: matrix building
# ---------------------------------------------------------------------------


class TestIntegrationMatrixBuilding:
    def test_matrix_shape(self, matrix, reader):
        n_taxa = len(reader.taxa)
        n_chars = len(reader.characters)
        assert matrix.shape == (n_taxa, n_chars)

    def test_matrix_index_matches_taxa_names(self, matrix, reader):
        expected_names = [t.name for t in reader.taxa]
        assert list(matrix.df.index) == expected_names

    def test_matrix_columns_match_character_ids(self, matrix, reader):
        expected_cols = [f"C{c.char_id}" for c in reader.characters]
        assert list(matrix.df.columns) == expected_cols

    def test_matrix_has_no_all_nan_row(self, matrix):
        """No taxon should be entirely missing."""
        all_nan_rows = matrix.df.index[matrix.df.isna().all(axis=1)]
        assert len(all_nan_rows) == 0


# ---------------------------------------------------------------------------
# Tests: matrix cleaning
# ---------------------------------------------------------------------------


class TestIntegrationMatrixCleaning:
    def test_remove_constant_does_not_increase_characters(self, matrix):
        from delta_phylo.matrix.matrix_cleaning import MatrixCleaner

        cleaned = MatrixCleaner(matrix).remove_constant_characters()
        assert cleaned.df.shape[1] <= matrix.df.shape[1]

    def test_filter_completeness_does_not_increase_shape(self, matrix):
        from delta_phylo.matrix.matrix_cleaning import MatrixCleaner

        cleaned = MatrixCleaner(matrix).filter_by_completeness(
            min_taxon_completeness=0.5, min_char_completeness=0.5
        )
        assert cleaned.df.shape[0] <= matrix.df.shape[0]
        assert cleaned.df.shape[1] <= matrix.df.shape[1]


# ---------------------------------------------------------------------------
# Tests: distance-based NJ tree
# ---------------------------------------------------------------------------


class TestIntegrationDistanceNJ:
    def test_nj_returns_tree(self, matrix):
        from delta_phylo.phylogeny.distance import DistanceAnalysis

        tree = DistanceAnalysis(matrix, metric="hamming").run()
        assert isinstance(tree, Tree)

    def test_nj_tree_has_correct_taxa(self, matrix):
        from delta_phylo.phylogeny.distance import DistanceAnalysis

        tree = DistanceAnalysis(matrix, metric="hamming").run()
        tree_taxa = {c.name for c in tree.get_terminals()}
        matrix_taxa = set(matrix.get_taxa_names())
        assert tree_taxa == matrix_taxa

    def test_simple_matching_nj_returns_tree(self, matrix):
        from delta_phylo.phylogeny.distance import DistanceAnalysis

        tree = DistanceAnalysis(matrix, metric="simple_matching").run()
        assert isinstance(tree, Tree)


# ---------------------------------------------------------------------------
# Tests: parsimony analysis
# ---------------------------------------------------------------------------


class TestIntegrationParsimony:
    def test_parsimony_returns_trees(self, matrix):
        from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis

        analysis = ParsimonyAnalysis(
            matrix=matrix, num_replicates=1, max_iterations=10, seed=42
        )
        trees = analysis.run()
        assert len(trees) >= 1

    def test_parsimony_best_score_non_negative(self, matrix):
        from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis

        analysis = ParsimonyAnalysis(
            matrix=matrix, num_replicates=1, max_iterations=10, seed=42
        )
        analysis.run()
        assert analysis.best_score >= 0

    def test_parsimony_tree_contains_all_taxa(self, matrix):
        from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis

        analysis = ParsimonyAnalysis(
            matrix=matrix, num_replicates=1, max_iterations=10, seed=42
        )
        trees = analysis.run()
        tree_taxa = {c.name for c in trees[0].get_terminals()}
        matrix_taxa = set(matrix.get_taxa_names())
        assert tree_taxa == matrix_taxa


# ---------------------------------------------------------------------------
# Tests: NEXUS export
# ---------------------------------------------------------------------------


class TestIntegrationNexusExport:
    def test_nexus_file_written(self, matrix, tmp_path):
        from delta_phylo.io.nexus_writer import NexusWriter

        out = tmp_path / "integration.nex"
        NexusWriter(matrix).write(out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_nexus_starts_with_header(self, matrix, tmp_path):
        from delta_phylo.io.nexus_writer import NexusWriter

        out = tmp_path / "integration.nex"
        NexusWriter(matrix).write(out)
        content = out.read_text()
        assert content.startswith("#NEXUS")

    def test_nexus_contains_all_taxa(self, matrix, tmp_path):
        from delta_phylo.io.nexus_writer import NexusWriter

        out = tmp_path / "integration.nex"
        NexusWriter(matrix).write(out)
        content = out.read_text()
        for name in matrix.get_taxa_names():
            safe = f"'{name}'" if " " in name else name
            assert safe in content, f"Taxon {name!r} not found in NEXUS output"

    def test_nexus_contains_missing_marker(self, matrix, tmp_path):
        from delta_phylo.io.nexus_writer import NexusWriter

        out = tmp_path / "integration.nex"
        NexusWriter(matrix).write(out)
        content = out.read_text()
        # The example dataset may or may not have missing data, but the format
        # declaration should always specify MISSING=?
        assert "MISSING=?" in content
