"""
Maximum Parsimony analysis using the Fitch algorithm.

Implements:
  - Fitch parsimony scoring of a tree
  - Nearest-Neighbor Interchange (NNI) heuristic tree search
  - Multiple-replicate searches with random starting trees
  - Consensus tree generation
"""

from __future__ import annotations

import logging
import random
from copy import deepcopy
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix
from delta_phylo.phylogeny.tree_utils import TreeUtils

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Fitch parsimony scoring
# ---------------------------------------------------------------------------


def _fitch_score_column(
    tree: Tree, states: Dict[str, Optional[int]]
) -> int:
    """Compute the Fitch parsimony score for a single character column.

    Args:
        tree: Rooted Bio.Phylo Tree.
        states: Mapping from terminal name to integer state (None = missing).

    Returns:
        Parsimony score (number of state changes).
    """
    score = 0
    clade_sets: Dict[int, set] = {}  # id(clade) → set of possible states

    def _postorder(clade: Clade) -> set:
        nonlocal score
        if clade.is_terminal():
            name = clade.name
            val = states.get(name)
            if val is None:
                # Missing: any state possible — use sentinel empty set
                return set()
            return {val}
        child_sets = []
        for child in clade.clades:
            cs = _postorder(child)
            child_sets.append(cs)

        # Remove empty sets (missing taxa) before intersection
        non_empty = [cs for cs in child_sets if cs]
        if not non_empty:
            return set()

        intersection = non_empty[0].copy()
        for cs in non_empty[1:]:
            intersection &= cs

        if intersection:
            return intersection
        # Union + increment score
        score += 1
        union = set()
        for cs in non_empty:
            union |= cs
        return union

    if tree.root:
        _postorder(tree.root)
    return score


def fitch_parsimony_score(tree: Tree, matrix: MorphologicalMatrix) -> int:
    """Compute the total Fitch parsimony score for a tree.

    Args:
        tree: A rooted Bio.Phylo Tree.
        matrix: MorphologicalMatrix providing character data.

    Returns:
        Total parsimony score (sum over all characters).
    """
    arr = matrix.to_numpy(missing_value=-1)
    taxa_names = matrix.get_taxa_names()
    total = 0

    for col_idx in range(arr.shape[1]):
        states: Dict[str, Optional[int]] = {}
        for row_idx, taxon_name in enumerate(taxa_names):
            val = arr[row_idx, col_idx]
            states[taxon_name] = None if val == -1 else int(val)
        total += _fitch_score_column(tree, states)

    return total


# ---------------------------------------------------------------------------
# NNI heuristic search
# ---------------------------------------------------------------------------


def _collect_internal_clades(tree: Tree) -> List[Clade]:
    """Return a list of internal (non-terminal) clades."""
    internals: List[Clade] = []

    def _walk(clade: Clade) -> None:
        if not clade.is_terminal():
            internals.append(clade)
            for child in clade.clades:
                _walk(child)

    _walk(tree.root)
    return internals


def _nni_neighbours(tree: Tree) -> List[Tree]:
    """Generate all NNI neighbours of *tree*.

    For each internal edge (parent → child), swap one subtree of the parent
    with one subtree of the child to produce two alternative trees.

    Args:
        tree: Source tree.

    Returns:
        List of new Tree objects (neighbours).
    """
    neighbours: List[Tree] = []
    internals = _collect_internal_clades(tree)

    for internal in internals:
        if len(internal.clades) < 2:
            continue
        parent = _find_parent(tree.root, internal)
        if parent is None:
            continue
        # Identify sibling subtrees of internal within parent
        siblings = [c for c in parent.clades if c is not internal]
        if not siblings:
            continue
        sibling = siblings[0]

        # NNI: swap each child of internal with sibling
        for i, child_of_internal in enumerate(internal.clades):
            new_tree = deepcopy(tree)
            # Re-locate nodes in the copy
            new_internal = _find_matching_clade(new_tree.root, internal)
            new_parent = _find_parent(new_tree.root, new_internal)
            if new_parent is None or new_internal is None:
                continue
            new_siblings = [c for c in new_parent.clades if c is not new_internal]
            if not new_siblings:
                continue
            new_sibling = new_siblings[0]
            new_child = new_internal.clades[i]

            # Perform swap
            idx_child = new_internal.clades.index(new_child)
            idx_sibling = new_parent.clades.index(new_sibling)
            new_internal.clades[idx_child] = new_sibling
            new_parent.clades[idx_sibling] = new_child
            neighbours.append(new_tree)

    return neighbours


def _find_parent(root: Clade, target: Clade) -> Optional[Clade]:
    """Find the parent clade of *target* in the tree rooted at *root*."""
    if target in root.clades:
        return root
    for child in root.clades:
        result = _find_parent(child, target)
        if result is not None:
            return result
    return None


def _find_matching_clade(root: Clade, original: Clade) -> Optional[Clade]:
    """Find a clade in a copied tree that structurally matches *original*."""
    # Match by terminal name sets
    orig_terminals = frozenset(c.name for c in original.get_terminals())

    def _search(clade: Clade) -> Optional[Clade]:
        terminals = frozenset(c.name for c in clade.get_terminals())
        if terminals == orig_terminals:
            return clade
        for child in clade.clades:
            result = _search(child)
            if result is not None:
                return result
        return None

    return _search(root)


def _random_tree(taxon_names: List[str], seed: Optional[int] = None) -> Tree:
    """Generate a random bifurcating tree from a list of taxon names.

    Uses a stepwise addition approach.

    Args:
        taxon_names: List of taxon names.
        seed: Optional random seed.

    Returns:
        Randomly generated Bio.Phylo Tree.
    """
    rng = random.Random(seed)
    names = list(taxon_names)
    rng.shuffle(names)

    if len(names) < 2:
        clades = [Clade(name=n) for n in names]
        return Tree(root=Clade(clades=clades), rooted=False)

    # Start with 3-taxon tree
    root = Clade(
        clades=[Clade(name=names[0]), Clade(name=names[1]), Clade(name=names[2])]
    )
    tree = Tree(root=root, rooted=False)

    # Add remaining taxa by random insertion
    for name in names[3:]:
        all_clades = list(tree.find_clades())
        attachment = rng.choice(all_clades)
        new_internal = Clade(clades=[Clade(name=name), deepcopy(attachment)])
        parent = _find_parent(tree.root, attachment)
        if parent is None:
            # attachment is root; wrap root
            new_root = Clade(clades=[Clade(name=name), tree.root])
            tree = Tree(root=new_root, rooted=False)
        else:
            idx = parent.clades.index(attachment)
            parent.clades[idx] = new_internal

    return tree


# ---------------------------------------------------------------------------
# ParsimonyAnalysis class
# ---------------------------------------------------------------------------


class ParsimonyAnalysis:
    """Maximum parsimony tree search using Fitch scoring and NNI.

    Args:
        matrix: MorphologicalMatrix to analyse.
        num_replicates: Number of random starting trees to use.
        max_iterations: Maximum NNI iterations per replicate.
        seed: Random seed for reproducibility.
    """

    def __init__(
        self,
        matrix: MorphologicalMatrix,
        num_replicates: int = 10,
        max_iterations: int = 100,
        seed: Optional[int] = 42,
    ) -> None:
        self.matrix = matrix
        self.num_replicates = num_replicates
        self.max_iterations = max_iterations
        self.seed = seed
        self.best_trees: List[Tree] = []
        self.best_score: int = int(1e9)

    def run(self) -> List[Tree]:
        """Execute the parsimony search.

        Returns:
            List of best-scoring trees found.
        """
        taxon_names = self.matrix.get_taxa_names()
        logger.info(
            "Starting parsimony search: %d replicates, %d max iter, %d taxa",
            self.num_replicates,
            self.max_iterations,
            len(taxon_names),
        )

        for rep in range(self.num_replicates):
            seed = None if self.seed is None else self.seed + rep
            start_tree = _random_tree(taxon_names, seed=seed)
            tree, score = self._nni_search(start_tree)
            logger.info("Replicate %d: score=%d", rep + 1, score)

            if score < self.best_score:
                self.best_score = score
                self.best_trees = [tree]
            elif score == self.best_score:
                self.best_trees.append(tree)

        logger.info(
            "Parsimony search complete. Best score=%d, trees found=%d",
            self.best_score,
            len(self.best_trees),
        )
        return self.best_trees

    def _nni_search(self, tree: Tree) -> Tuple[Tree, int]:
        """NNI heuristic search from a single starting tree.

        Args:
            tree: Starting tree.

        Returns:
            Tuple of (best tree, best score).
        """
        current_tree = deepcopy(tree)
        current_score = fitch_parsimony_score(current_tree, self.matrix)

        for iteration in range(self.max_iterations):
            improved = False
            neighbours = _nni_neighbours(current_tree)
            for neighbour in neighbours:
                score = fitch_parsimony_score(neighbour, self.matrix)
                if score < current_score:
                    current_tree = neighbour
                    current_score = score
                    improved = True
                    break  # restart from improved tree
            if not improved:
                break  # local optimum reached

        return current_tree, current_score

    def consensus_tree(self, cutoff: float = 0.5) -> Tree:
        """Return a majority-rule consensus of all best trees.

        Args:
            cutoff: Minimum frequency for clade inclusion (0–1).

        Returns:
            Consensus Bio.Phylo Tree.
        """
        if not self.best_trees:
            raise RuntimeError("No trees found. Run .run() first.")
        return TreeUtils.majority_rule_consensus(self.best_trees, cutoff=cutoff)
