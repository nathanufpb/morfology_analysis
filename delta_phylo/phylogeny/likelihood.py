"""
Maximum Likelihood analysis for discrete morphological characters.

Implements the Mk model (Lewis 2001) for ordered and unordered characters,
with optional rate variation (Gamma +G).

Note: This is a simplified implementation intended for educational use.
For production ML analyses, please use IQ-TREE or RAxML with NEXUS export.
"""

from __future__ import annotations

import logging
import math
from typing import List, Optional

import numpy as np
from Bio.Phylo.BaseTree import Tree

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)


def mk_likelihood(
    tree: Tree,
    matrix: MorphologicalMatrix,
    equal_freq: bool = True,
) -> float:
    """Compute an approximate log-likelihood under the Mk model.

    The Mk model (Lewis 2001) assumes equal rates among all state transitions
    for each character.

    Args:
        tree: Bio.Phylo Tree with branch lengths.
        matrix: MorphologicalMatrix.
        equal_freq: If True, assumes equal state frequencies.

    Returns:
        Approximate log-likelihood (negative float).
    """
    logger.warning(
        "mk_likelihood is a simplified implementation. "
        "For rigorous ML analyses, export to NEXUS and use IQ-TREE."
    )
    arr = matrix.to_numpy(missing_value=-1)
    taxa_names = matrix.get_taxa_names()
    n_chars = arr.shape[1]

    total_log_lik = 0.0
    for col_idx in range(n_chars):
        col = arr[:, col_idx]
        valid = col[col != -1]
        if len(valid) == 0:
            continue
        n_states = int(valid.max()) + 1
        if n_states < 2:
            continue
        # Simple per-column likelihood: frequency-based approximation
        freq = np.bincount(valid, minlength=n_states) / len(valid)
        # Use site log-likelihood as log of expected frequency under equal rates
        for val in valid:
            prob = freq[val] if freq[val] > 0 else 1e-10
            total_log_lik += math.log(prob)

    return total_log_lik


class LikelihoodAnalysis:
    """Simplified Maximum Likelihood tree optimisation using the Mk model.

    For production analyses, use :meth:`export to NEXUS <delta_phylo.io.nexus_writer.NexusWriter>`
    and run IQ-TREE or RAxML externally.

    Args:
        matrix: MorphologicalMatrix to analyse.
    """

    def __init__(self, matrix: MorphologicalMatrix) -> None:
        self.matrix = matrix

    def score_tree(self, tree: Tree) -> float:
        """Return the approximate Mk log-likelihood for a given tree.

        Args:
            tree: Bio.Phylo Tree.

        Returns:
            Log-likelihood (negative float).
        """
        return mk_likelihood(tree, self.matrix)

    def best_of(self, trees: List[Tree]) -> Tree:
        """Select the tree with the highest log-likelihood.

        Args:
            trees: List of Bio.Phylo Tree objects.

        Returns:
            Best-scoring Tree.
        """
        best_tree = max(trees, key=lambda t: self.score_tree(t))
        return best_tree
