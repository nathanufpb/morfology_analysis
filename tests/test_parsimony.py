"""Tests for parsimony analysis."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
from Bio.Phylo.BaseTree import Tree

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder, MorphologicalMatrix
from delta_phylo.phylogeny.parsimony import (
    ParsimonyAnalysis,
    fitch_parsimony_score,
    _random_tree,
    _nni_neighbours,
)
from delta_phylo.phylogeny.distance import DistanceAnalysis, hamming_distance, simple_matching
from delta_phylo.phylogeny.tree_utils import TreeUtils


# ---------------------------------------------------------------------------
# Fixture: small morphological matrix
# ---------------------------------------------------------------------------


@pytest.fixture
def small_matrix() -> MorphologicalMatrix:
    """Build a 5-taxon × 6-character matrix."""
    chars = textwrap.dedent("""\
        #1. A/ 1. 0/ 2. 1/
        #2. B/ 1. 0/ 2. 1/
        #3. C/ 1. 0/ 2. 1/ 3. 2/
        #4. D/ 1. 0/ 2. 1/
        #5. E/ 1. 0/ 2. 1/
        #6. F/ 1. 0/ 2. 1/
    """)
    items = textwrap.dedent("""\
        #1. Alpha/ 1,1 2,1 3,1 4,1 5,1 6,1/
        #2. Beta/  1,2 2,2 3,2 4,1 5,1 6,1/
        #3. Gamma/ 1,2 2,2 3,3 4,2 5,2 6,1/
        #4. Delta/ 1,1 2,1 3,2 4,2 5,2 6,2/
        #5. Epsilon/ 1,1 2,2 3,1 4,1 5,2 6,2/
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
# Tests: Fitch parsimony score
# ---------------------------------------------------------------------------


class TestFitchParsimony:
    def test_score_non_negative(self, small_matrix: MorphologicalMatrix) -> None:
        tree = _random_tree(small_matrix.get_taxa_names(), seed=0)
        # Ensure tree is rooted for Fitch algorithm
        tree.rooted = True
        score = fitch_parsimony_score(tree, small_matrix)
        assert score >= 0

    def test_score_type_is_int(self, small_matrix: MorphologicalMatrix) -> None:
        tree = _random_tree(small_matrix.get_taxa_names(), seed=1)
        score = fitch_parsimony_score(tree, small_matrix)
        assert isinstance(score, int)


# ---------------------------------------------------------------------------
# Tests: ParsimonyAnalysis
# ---------------------------------------------------------------------------


class TestParsimonyAnalysis:
    def test_run_returns_trees(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = ParsimonyAnalysis(
            matrix=small_matrix, num_replicates=3, max_iterations=20, seed=42
        )
        trees = analysis.run()
        assert len(trees) >= 1
        assert all(isinstance(t, Tree) for t in trees)

    def test_best_score_non_negative(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = ParsimonyAnalysis(
            matrix=small_matrix, num_replicates=3, max_iterations=20, seed=42
        )
        analysis.run()
        assert analysis.best_score >= 0

    def test_all_taxa_in_tree(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = ParsimonyAnalysis(
            matrix=small_matrix, num_replicates=2, max_iterations=10, seed=0
        )
        trees = analysis.run()
        expected_taxa = set(small_matrix.get_taxa_names())
        for tree in trees:
            found = set(TreeUtils.terminal_names(tree))
            assert found == expected_taxa

    def test_consensus_tree(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = ParsimonyAnalysis(
            matrix=small_matrix, num_replicates=5, max_iterations=20, seed=7
        )
        analysis.run()
        consensus = analysis.consensus_tree(cutoff=0.5)
        assert consensus is not None


# ---------------------------------------------------------------------------
# Tests: Distance analysis
# ---------------------------------------------------------------------------


class TestDistanceFunctions:
    def test_hamming_identical(self) -> None:
        import numpy as np

        a = np.array([0, 1, 2, 1, 0])
        assert hamming_distance(a, a) == 0.0

    def test_hamming_completely_different(self) -> None:
        import numpy as np

        a = np.array([0, 0, 0, 0])
        b = np.array([1, 1, 1, 1])
        assert hamming_distance(a, b) == 1.0

    def test_hamming_with_missing(self) -> None:
        import numpy as np

        a = np.array([0, -1, 1, -1])
        b = np.array([0, 1, 0, -1])
        # Only positions 0 and 2 are comparable: 0==0 (match), 1!=0 (mismatch)
        d = hamming_distance(a, b)
        assert d == pytest.approx(0.5)

    def test_simple_matching_identical(self) -> None:
        import numpy as np

        a = np.array([0, 1, 2, 1])
        assert simple_matching(a, a) == 0.0

    def test_hamming_all_missing(self) -> None:
        import numpy as np

        a = np.array([-1, -1])
        b = np.array([-1, -1])
        assert hamming_distance(a, b) == 1.0


class TestDistanceAnalysis:
    def test_nj_returns_tree(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = DistanceAnalysis(matrix=small_matrix, metric="hamming")
        tree = analysis.run()
        assert isinstance(tree, Tree)

    def test_nj_tree_has_all_taxa(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = DistanceAnalysis(matrix=small_matrix)
        tree = analysis.run()
        expected = set(small_matrix.get_taxa_names())
        found = set(TreeUtils.terminal_names(tree))
        assert found == expected

    def test_nj_simple_matching(self, small_matrix: MorphologicalMatrix) -> None:
        analysis = DistanceAnalysis(matrix=small_matrix, metric="simple_matching")
        tree = analysis.run()
        assert isinstance(tree, Tree)


# ---------------------------------------------------------------------------
# Tests: TreeUtils
# ---------------------------------------------------------------------------


class TestTreeUtils:
    def test_random_tree_taxa_count(self) -> None:
        names = ["A", "B", "C", "D", "E"]
        tree = _random_tree(names, seed=0)
        found = TreeUtils.terminal_names(tree)
        assert set(found) == set(names)

    def test_newick_roundtrip(self) -> None:
        names = ["Alpha", "Beta", "Gamma", "Delta"]
        tree = _random_tree(names, seed=42)
        newick = TreeUtils.tree_to_newick(tree)
        tree2 = TreeUtils.newick_to_tree(newick)
        found = set(TreeUtils.terminal_names(tree2))
        assert found == set(names)

    def test_star_tree(self) -> None:
        names = ["T1", "T2", "T3"]
        tree = TreeUtils.star_tree(names)
        assert set(TreeUtils.terminal_names(tree)) == set(names)

    def test_save_load_tree(self, tmp_path: Path) -> None:
        names = ["P", "Q", "R", "S"]
        tree = _random_tree(names, seed=10)
        out = tmp_path / "test.nwk"
        TreeUtils.save_tree(tree, str(out))
        loaded = TreeUtils.load_tree(str(out))
        assert set(TreeUtils.terminal_names(loaded)) == set(names)
