"""
Maximum Likelihood analysis for discrete morphological characters.

Implements the Mk model (Lewis 2001) for ordered and unordered characters,
using the Felsenstein pruning algorithm to compute site log-likelihoods that
depend on the tree topology and branch lengths.

Note: This is a simplified implementation intended for educational use.
For production ML analyses, please use IQ-TREE or RAxML with NEXUS export.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix

logger = logging.getLogger(__name__)

# Default branch length used when a clade has none specified.
_DEFAULT_BRANCH_LENGTH = 0.1


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _mk_transition_matrix(n_states: int, branch_length: float, mu: float = 1.0) -> np.ndarray:
    """Compute the Mk-model transition probability matrix P(t) = exp(Q*t).

    For the Mk model all off-diagonal rates are equal to *mu* and the
    diagonal entries q_ii = -(n_states - 1) * mu.  The closed-form solution is

        P(t) = J/k + (I - J/k) * exp(-k * mu * t)

    where J is the k×k all-ones matrix, I is the identity and k = n_states.

    Args:
        n_states: Number of character states k.
        branch_length: Edge length t (expected substitutions per site).
        mu: Off-diagonal rate parameter (default 1.0).

    Returns:
        k×k numpy array of transition probabilities.
    """
    k = n_states
    t = max(branch_length, 1e-9)   # guard against zero / negative lengths
    e_neg = math.exp(-k * mu * t)
    ones_over_k = np.full((k, k), 1.0 / k)
    P = ones_over_k + (np.eye(k) - ones_over_k) * e_neg
    return P


def _pruning_partial(
    clade: Clade,
    state_lookup: Dict[str, int],
    n_states: int,
    mu: float = 1.0,
) -> np.ndarray:
    """Recursively compute Felsenstein partial likelihoods for *clade*.

    Args:
        clade: Current Bio.Phylo Clade (may be leaf or internal).
        state_lookup: Mapping from taxon name to observed state index (-1 = missing).
        n_states: Number of character states k.
        mu: Mk-model rate parameter.

    Returns:
        1-D numpy array of shape (n_states,) with partial likelihoods.
    """
    if clade.is_terminal():
        partial = np.zeros(n_states)
        state = state_lookup.get(clade.name, -1)
        if state == -1:           # missing data → uniform over all states
            partial[:] = 1.0
        else:
            partial[int(state)] = 1.0
        return partial

    # Internal node: product over children of summed transition contributions.
    partial = np.ones(n_states)
    for child in clade.clades:
        bl = child.branch_length if child.branch_length is not None else _DEFAULT_BRANCH_LENGTH
        P = _mk_transition_matrix(n_states, bl, mu)           # (k, k)
        child_partial = _pruning_partial(child, state_lookup, n_states, mu)
        # For each parent state s: contribution = Σ_t P[s,t] * L_child[t]
        partial *= P @ child_partial
    return partial


def mk_likelihood(
    tree: Tree,
    matrix: MorphologicalMatrix,
    equal_freq: bool = True,
    mu: float = 1.0,
) -> float:
    """Compute log-likelihood under the Mk model via Felsenstein pruning.

    Traverses the tree in a post-order fashion (pruning algorithm, Felsenstein
    1981) to compute the exact site likelihood for each character column.
    The score therefore depends on the tree topology and branch lengths.

    Args:
        tree: Bio.Phylo Tree with branch lengths.
        matrix: MorphologicalMatrix.
        equal_freq: If True (default) uses equal state-frequency prior π_i = 1/k.
        mu: Off-diagonal substitution rate for the Mk model.

    Returns:
        Log-likelihood (a negative float; higher = better fit).
    """
    arr = matrix.to_numpy(missing_value=-1)
    taxa_names = matrix.get_taxa_names()
    n_chars = arr.shape[1]

    # Build taxon-index → name map used inside the per-character loop.
    idx_to_name: List[str] = list(taxa_names)

    total_log_lik = 0.0

    for col_idx in range(n_chars):
        col = arr[:, col_idx]
        valid = col[col != -1]
        if len(valid) == 0:
            continue
        n_states = int(valid.max()) + 1
        if n_states < 2:
            continue

        # Build per-taxon state lookup for this character column.
        state_lookup: Dict[str, int] = {
            idx_to_name[i]: int(col[i])
            for i in range(len(col))
            if col[i] != -1
        }
        # missing taxa get -1 (handled inside _pruning_partial as uniform)

        # Felsenstein pruning from the root.
        root_partial = _pruning_partial(tree.root, state_lookup, n_states, mu)

        # Site likelihood = Σ_s π(s) * L_root(s)
        prior = 1.0 / n_states   # equal frequency
        site_lik = float(prior * root_partial.sum())
        if site_lik <= 0.0:
            site_lik = 1e-300     # numerical floor to avoid log(0)

        total_log_lik += math.log(site_lik)

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
