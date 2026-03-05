"""
Distance-based phylogenetic analysis.

Implements:
  - Hamming distance
  - Simple matching coefficient
  - Gower distance (handles numeric/continuous characters)
  - Neighbor Joining tree inference using Bio.Phylo
"""

from __future__ import annotations

import logging
from typing import List, Optional

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from delta_phylo.matrix.matrix_builder import MorphologicalMatrix
from delta_phylo.parser.characters import CharacterType
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


def gower_distance(
    row_a: np.ndarray,
    row_b: np.ndarray,
    is_numeric: np.ndarray,
    ranges: np.ndarray,
) -> float:
    """Compute Gower's distance between two character rows.

    Gower's coefficient handles mixed data (numeric + discrete).  For each
    character *j*:

    * **Numeric** (continuous): contribution = |a_j − b_j| / range_j
    * **Discrete**: contribution = 0 if a_j == b_j, else 1

    Missing values (NaN for numeric, −1 for discrete) are excluded.

    Args:
        row_a: Float array for taxon A (NaN = missing for numeric chars).
        row_b: Float array for taxon B.
        is_numeric: Boolean array indicating which columns are numeric.
        ranges: For numeric columns: the pre-computed per-column value range
            (max − min).  Set to 1.0 for discrete columns / zero-range columns.

    Returns:
        Gower distance in [0, 1].  Returns 1.0 if no comparable characters.
    """
    total = 0.0
    diffs = 0.0
    n = len(row_a)
    for j in range(n):
        a_j = row_a[j]
        b_j = row_b[j]
        if is_numeric[j]:
            if np.isnan(a_j) or np.isnan(b_j):
                continue
            r = ranges[j] if ranges[j] > 0 else 1.0
            diffs += abs(a_j - b_j) / r
        else:
            if int(a_j) == -1 or int(b_j) == -1:
                continue
            diffs += 0.0 if a_j == b_j else 1.0
        total += 1.0
    if total == 0:
        return 1.0
    return diffs / total


def _build_distance_matrix(
    matrix: MorphologicalMatrix, metric: str = "hamming"
) -> np.ndarray:
    """Build a pairwise distance matrix from a morphological matrix.

    For matrices containing numeric characters, the ``"gower"`` metric is
    automatically used regardless of the *metric* argument.

    Args:
        matrix: MorphologicalMatrix object.
        metric: ``"hamming"``, ``"simple_matching"``, or ``"gower"``.

    Returns:
        Square NumPy distance matrix of shape (n_taxa, n_taxa).
    """
    # Detect if any character is numeric
    has_numeric = any(
        c.char_type == CharacterType.NUMERIC for c in matrix.characters
    )

    if metric == "gower" or has_numeric:
        return _build_gower_matrix(matrix)

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


def _build_gower_matrix(matrix: MorphologicalMatrix) -> np.ndarray:
    """Build a pairwise Gower-distance matrix.

    Numeric columns keep their raw float values; discrete columns are encoded
    as integers with −1 for missing.

    Args:
        matrix: MorphologicalMatrix.

    Returns:
        Square float NumPy distance array.
    """
    float_arr = matrix.to_numpy_float(missing_value=float("nan"))
    int_arr = matrix.to_numpy(missing_value=-1)

    is_numeric = np.array(
        [c.char_type == CharacterType.NUMERIC for c in matrix.characters],
        dtype=bool,
    )

    # Build the combined array: numeric columns → float values, others → int-encoded -1 for missing
    combined = float_arr.copy()
    for j in range(combined.shape[1]):
        if not is_numeric[j]:
            # Use the int-encoded column (handles NaN → -1)
            combined[:, j] = int_arr[:, j].astype(float)

    # Compute per-column ranges for numeric columns
    ranges = np.ones(combined.shape[1])
    for j in range(combined.shape[1]):
        if is_numeric[j]:
            col = combined[:, j]
            col_valid = col[~np.isnan(col)]
            if len(col_valid) > 1:
                ranges[j] = float(col_valid.max() - col_valid.min())

    n = combined.shape[0]
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = gower_distance(combined[i], combined[j], is_numeric, ranges)
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
