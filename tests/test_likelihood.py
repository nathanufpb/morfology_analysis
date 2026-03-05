"""Tests for LikelihoodAnalysis and mk_likelihood — closes issue #9.

NOTE: The current mk_likelihood implementation intentionally ignores the tree
topology (see issue #4). Tests documenting this known-broken behaviour are
marked with `xfail` so they serve as regression tests after the fix.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest
from Bio.Phylo.BaseTree import Tree

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder, MorphologicalMatrix
from delta_phylo.phylogeny.likelihood import mk_likelihood, LikelihoodAnalysis
from delta_phylo.phylogeny.tree_utils import TreeUtils
from delta_phylo.phylogeny.parsimony import _random_tree


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_matrix() -> MorphologicalMatrix:
    chars = textwrap.dedent("""\
        #1. A/ 1. 0/ 2. 1/
        #2. B/ 1. 0/ 2. 1/
        #3. C/ 1. 0/ 2. 1/
        #4. D/ 1. 0/ 2. 1/
    """)
    items = textwrap.dedent("""\
        #1. Alpha/   1,1 2,1 3,1 4,1/
        #2. Beta/    1,2 2,2 3,1 4,1/
        #3. Gamma/   1,2 2,2 3,2 4,1/
        #4. Delta/   1,1 2,1 3,2 4,2/
    """)
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp)
        (p / "characters").write_text(chars)
        (p / "items").write_text(items)
        reader = DeltaReader(p)
        reader.read()
        return MatrixBuilder(reader.characters, reader.taxa).build()


@pytest.fixture
def tree_a(small_matrix) -> Tree:
    """A random tree over the small_matrix taxa."""
    return _random_tree(small_matrix.get_taxa_names(), seed=0)


@pytest.fixture
def tree_b(small_matrix) -> Tree:
    """A different random tree over the same taxa."""
    return _random_tree(small_matrix.get_taxa_names(), seed=99)


# ---------------------------------------------------------------------------
# Tests: mk_likelihood (functional properties)
# ---------------------------------------------------------------------------


class TestMkLikelihood:
    def test_returns_negative_float(self, tree_a, small_matrix):
        ll = mk_likelihood(tree_a, small_matrix)
        assert isinstance(ll, float)
        assert ll <= 0.0

    def test_returns_finite_value(self, tree_a, small_matrix):
        ll = mk_likelihood(tree_a, small_matrix)
        assert math.isfinite(ll)

    def test_equal_freq_flag_accepted(self, tree_a, small_matrix):
        """equal_freq=True and False should both return valid floats."""
        ll_true = mk_likelihood(tree_a, small_matrix, equal_freq=True)
        ll_false = mk_likelihood(tree_a, small_matrix, equal_freq=False)
        assert math.isfinite(ll_true)
        assert math.isfinite(ll_false)

    @pytest.mark.xfail(
        reason="Issue #4: mk_likelihood ignores tree topology; "
               "expected to fail until Felsenstein pruning is implemented."
    )
    def test_different_trees_get_different_scores(self, tree_a, tree_b, small_matrix):
        """After issue #4 is fixed, different trees must receive different scores."""
        ll_a = mk_likelihood(tree_a, small_matrix)
        ll_b = mk_likelihood(tree_b, small_matrix)
        assert ll_a != ll_b


# ---------------------------------------------------------------------------
# Tests: LikelihoodAnalysis
# ---------------------------------------------------------------------------


class TestLikelihoodAnalysis:
    def test_score_tree_returns_float(self, tree_a, small_matrix):
        analysis = LikelihoodAnalysis(small_matrix)
        score = analysis.score_tree(tree_a)
        assert isinstance(score, float)

    def test_best_of_returns_a_tree(self, tree_a, tree_b, small_matrix):
        analysis = LikelihoodAnalysis(small_matrix)
        best = analysis.best_of([tree_a, tree_b])
        assert isinstance(best, Tree)

    def test_best_of_single_tree(self, tree_a, small_matrix):
        analysis = LikelihoodAnalysis(small_matrix)
        best = analysis.best_of([tree_a])
        assert best is tree_a

    def test_best_of_empty_list_raises(self, small_matrix):
        analysis = LikelihoodAnalysis(small_matrix)
        with pytest.raises((ValueError, IndexError)):
            analysis.best_of([])

    @pytest.mark.xfail(
        reason="Issue #4: all trees score identically until topology is used."
    )
    def test_best_of_selects_higher_scoring_tree(self, tree_a, tree_b, small_matrix):
        """After issue #4 is fixed, best_of must pick the truly best tree."""
        analysis = LikelihoodAnalysis(small_matrix)
        ll_a = analysis.score_tree(tree_a)
        ll_b = analysis.score_tree(tree_b)
        best = analysis.best_of([tree_a, tree_b])
        expected_best = tree_a if ll_a > ll_b else tree_b
        # Compare Newick strings as a proxy for topology equality
        assert TreeUtils.tree_to_newick(best) == TreeUtils.tree_to_newick(expected_best)
