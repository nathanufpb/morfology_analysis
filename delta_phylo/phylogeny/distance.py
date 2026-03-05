"""
Distance-based phylogenetic analysis.

Implements:
  - Hamming distance
  - Simple matching coefficient
  - Neighbor Joining tree inference using Bio.Phylo
"""

from __future__ import annotations

import logging
from typing import List

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix
from delta_phylo.phylogeny.tree_utils import TreeUtils

logger = logging.getLogger(__name__)


def hamming_distance(row_a: np.ndarray, row_b: np.ndarray) -> float:
    """Compute the normalised Hamming distance between two character rows.

    Missing values (encoded as ``-1``) are excluded from the calculation.

    Args:
        row_a: Integer array for taxon A.
        row_b: Integer array for taxon B.

    Returns:
        Normalised Hamming distance in [0, 1].  Returns 1.0 if no
        comparable characters exist.
    """
    mask = (row_a != -1) & (row_b != -1)
    if mask.sum() == 0:
        return 1.0
    return float((row_a[mask] != row_b[mask]).sum()) / float(mask.sum())


def simple_matching(row_a: np.ndarray, row_b: np.ndarray) -> float:
    """Compute the simple matching dissimilarity between two character rows.

    Missing values (encoded as ``-1``) are excluded.

    Args:
        row_a: Integer array for taxon A.
        row_b: Integer array for taxon B.

    Returns:
        Proportion of characters *not* in common (1 − SMC).
    """
    mask = (row_a != -1) & (row_b != -1)
    if mask.sum() == 0:
        return 1.0
    matches = float((row_a[mask] == row_b[mask]).sum())
    return 1.0 - matches / float(mask.sum())


def _build_distance_matrix(
    matrix: MorphologicalMatrix, metric: str = "hamming"
) -> np.ndarray:
    """Build a pairwise distance matrix from a morphological matrix.

    Args:
        matrix: MorphologicalMatrix object.
        metric: ``"hamming"`` or ``"simple_matching"``.

    Returns:
        Square NumPy distance matrix of shape (n_taxa, n_taxa).
    """
    arr = matrix.to_numpy(missing_value=-1)
    n = arr.shape[0]
    dist = np.zeros((n, n))

    if metric == "hamming":
        fn = hamming_distance
    elif metric == "simple_matching":
        fn = simple_matching
    else:
        raise ValueError(f"Unknown metric: {metric!r}")

    for i in range(n):
        for j in range(i + 1, n):
            d = fn(arr[i], arr[j])
            dist[i, j] = d
            dist[j, i] = d

    return dist


class DistanceAnalysis:
    """Distance-based tree inference using Neighbor Joining.

    Args:
        matrix: MorphologicalMatrix to analyse.
        metric: Distance metric to use (``"hamming"`` or ``"simple_matching"``).
    """

    def __init__(
        self,
        matrix: MorphologicalMatrix,
        metric: str = "hamming",
    ) -> None:
        self.matrix = matrix
        self.metric = metric

    def run(self):
        """Infer a Neighbor Joining tree.

        Returns:
            Bio.Phylo Tree object.
        """
        taxa_names = self.matrix.get_taxa_names()
        logger.info(
            "Running Neighbor Joining with metric=%r on %d taxa",
            self.metric,
            len(taxa_names),
        )
        dist_arr = _build_distance_matrix(self.matrix, metric=self.metric)
        biopython_dm = self._to_biopython_dm(dist_arr, taxa_names)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(biopython_dm)
        logger.info("NJ tree inference complete.")
        return tree

    @staticmethod
    def _to_biopython_dm(
        dist_arr: np.ndarray, names: List[str]
    ) -> DistanceMatrix:
        """Convert a NumPy distance matrix to a Bio.Phylo DistanceMatrix.

        Bio.Phylo's DistanceMatrix requires a lower-triangular list-of-lists.

        Args:
            dist_arr: Square distance matrix.
            names: Taxon names.

        Returns:
            Bio.Phylo DistanceMatrix.
        """
        n = len(names)
        matrix_list = []
        for i in range(n):
            row = [dist_arr[i, j] for j in range(i + 1)]
            matrix_list.append(row)
        return DistanceMatrix(names=names, matrix=matrix_list)
