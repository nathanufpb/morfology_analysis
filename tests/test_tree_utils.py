"""Tests for TreeUtils — closes issue #11.

Covers:
  - newick_to_tree()
  - tree_to_newick()
  - round-trip Newick
  - save_tree / load_tree
  - majority_rule_consensus()
  - root_midpoint()
  - terminal_names()
  - star_tree()
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Phylo.BaseTree import Tree, Clade

from delta_phylo.phylogeny.tree_utils import TreeUtils


# ---------------------------------------------------------------------------
# Fixtures / shared Newick strings
# ---------------------------------------------------------------------------


SIMPLE_NEWICK = "((Alpha:0.1,Beta:0.2):0.1,(Gamma:0.3,Delta:0.4):0.2);"
TAXA = ["Alpha", "Beta", "Gamma", "Delta"]


@pytest.fixture
def simple_tree() -> Tree:
    return TreeUtils.newick_to_tree(SIMPLE_NEWICK)


# ---------------------------------------------------------------------------
# Tests: newick_to_tree
# ---------------------------------------------------------------------------


class TestNewickToTree:
    def test_returns_tree_object(self):
        tree = TreeUtils.newick_to_tree(SIMPLE_NEWICK)
        assert isinstance(tree, Tree)

    def test_correct_taxon_count(self):
        tree = TreeUtils.newick_to_tree(SIMPLE_NEWICK)
        assert len(tree.get_terminals()) == 4

    def test_single_leaf_newick(self):
        tree = TreeUtils.newick_to_tree("Alpha;")
        assert len(tree.get_terminals()) == 1

    def test_taxon_names_present(self):
        tree = TreeUtils.newick_to_tree(SIMPLE_NEWICK)
        names = {c.name for c in tree.get_terminals()}
        assert names == {"Alpha", "Beta", "Gamma", "Delta"}


# ---------------------------------------------------------------------------
# Tests: tree_to_newick
# ---------------------------------------------------------------------------


class TestTreeToNewick:
    def test_returns_string(self, simple_tree):
        result = TreeUtils.tree_to_newick(simple_tree)
        assert isinstance(result, str)

    def test_nonempty_string(self, simple_tree):
        result = TreeUtils.tree_to_newick(simple_tree)
        assert len(result) > 0

    def test_round_trip_preserves_taxa(self):
        tree1 = TreeUtils.newick_to_tree(SIMPLE_NEWICK)
        nwk = TreeUtils.tree_to_newick(tree1)
        tree2 = TreeUtils.newick_to_tree(nwk)
        names1 = {c.name for c in tree1.get_terminals()}
        names2 = {c.name for c in tree2.get_terminals()}
        assert names1 == names2

    def test_round_trip_preserves_topology(self):
        """Re-parsing the Newick should yield the same number of internal nodes."""
        tree1 = TreeUtils.newick_to_tree(SIMPLE_NEWICK)
        nwk = TreeUtils.tree_to_newick(tree1)
        tree2 = TreeUtils.newick_to_tree(nwk)
        assert len(list(tree1.get_nonterminals())) == len(
            list(tree2.get_nonterminals())
        )


# ---------------------------------------------------------------------------
# Tests: save_tree / load_tree
# ---------------------------------------------------------------------------


class TestSaveLoadTree:
    def test_save_creates_file(self, simple_tree, tmp_path):
        outfile = str(tmp_path / "tree.nwk")
        TreeUtils.save_tree(simple_tree, outfile, fmt="newick")
        assert Path(outfile).exists()
        assert Path(outfile).stat().st_size > 0

    def test_load_returns_tree(self, simple_tree, tmp_path):
        outfile = str(tmp_path / "tree.nwk")
        TreeUtils.save_tree(simple_tree, outfile, fmt="newick")
        loaded = TreeUtils.load_tree(outfile, fmt="newick")
        assert isinstance(loaded, Tree)

    def test_save_load_round_trip_taxa(self, simple_tree, tmp_path):
        outfile = str(tmp_path / "tree.nwk")
        TreeUtils.save_tree(simple_tree, outfile)
        loaded = TreeUtils.load_tree(outfile)
        original_names = {c.name for c in simple_tree.get_terminals()}
        loaded_names = {c.name for c in loaded.get_terminals()}
        assert original_names == loaded_names


# ---------------------------------------------------------------------------
# Tests: majority_rule_consensus
# ---------------------------------------------------------------------------


class TestMajorityRuleConsensus:
    def test_consensus_of_identical_trees_is_a_tree(self, simple_tree):
        trees = [TreeUtils.newick_to_tree(SIMPLE_NEWICK) for _ in range(5)]
        consensus = TreeUtils.majority_rule_consensus(trees, cutoff=0.5)
        assert isinstance(consensus, Tree)

    def test_consensus_contains_all_terminals(self, simple_tree):
        trees = [TreeUtils.newick_to_tree(SIMPLE_NEWICK) for _ in range(3)]
        consensus = TreeUtils.majority_rule_consensus(trees, cutoff=0.5)
        names = {c.name for c in consensus.get_terminals()}
        assert names == {"Alpha", "Beta", "Gamma", "Delta"}


# ---------------------------------------------------------------------------
# Tests: terminal_names
# ---------------------------------------------------------------------------


class TestTerminalNames:
    def test_returns_all_taxa(self, simple_tree):
        names = TreeUtils.terminal_names(simple_tree)
        assert sorted(names) == sorted(TAXA)

    def test_returns_list(self, simple_tree):
        assert isinstance(TreeUtils.terminal_names(simple_tree), list)


# ---------------------------------------------------------------------------
# Tests: star_tree
# ---------------------------------------------------------------------------


class TestStarTree:
    def test_returns_tree(self):
        tree = TreeUtils.star_tree(TAXA)
        assert isinstance(tree, Tree)

    def test_all_taxa_are_terminals(self):
        tree = TreeUtils.star_tree(TAXA)
        names = {c.name for c in tree.get_terminals()}
        assert names == set(TAXA)

    def test_single_internal_node(self):
        tree = TreeUtils.star_tree(TAXA)
        # Star tree has exactly one internal node (the root)
        non_terminals = list(tree.get_nonterminals())
        assert len(non_terminals) == 1


# ---------------------------------------------------------------------------
# Tests: root_midpoint
# ---------------------------------------------------------------------------


class TestRootMidpoint:
    def test_returns_tree(self, simple_tree):
        result = TreeUtils.root_midpoint(simple_tree)
        assert isinstance(result, Tree)
