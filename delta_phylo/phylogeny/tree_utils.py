"""
Tree utility functions using Bio.Phylo from Biopython.
"""

from __future__ import annotations

import io
import logging
from typing import List, Optional

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

logger = logging.getLogger(__name__)


class TreeUtils:
    """Utility class for phylogenetic tree operations.

    Provides static methods for tree manipulation, Newick I/O, and
    consensus tree computation.
    """

    @staticmethod
    def newick_to_tree(newick: str) -> Tree:
        """Parse a Newick string into a Bio.Phylo Tree object.

        Args:
            newick: Newick format string.

        Returns:
            Bio.Phylo Tree object.
        """
        handle = io.StringIO(newick)
        return Phylo.read(handle, "newick")

    @staticmethod
    def tree_to_newick(tree: Tree) -> str:
        """Convert a Bio.Phylo Tree to a Newick string.

        Args:
            tree: Bio.Phylo Tree object.

        Returns:
            Newick string representation.
        """
        handle = io.StringIO()
        Phylo.write(tree, handle, "newick")
        return handle.getvalue().strip()

    @staticmethod
    def save_tree(tree: Tree, filepath: str, fmt: str = "newick") -> None:
        """Write a tree to a file.

        Args:
            tree: Bio.Phylo Tree object.
            filepath: Path to output file.
            fmt: Format string accepted by Bio.Phylo (``"newick"``, ``"nexus"``).
        """
        Phylo.write(tree, filepath, fmt)
        logger.info("Tree saved to %s (format=%s)", filepath, fmt)

    @staticmethod
    def load_tree(filepath: str, fmt: str = "newick") -> Tree:
        """Load a tree from a file.

        Args:
            filepath: Path to input file.
            fmt: Format string accepted by Bio.Phylo.

        Returns:
            Bio.Phylo Tree object.
        """
        tree = Phylo.read(filepath, fmt)
        logger.info("Tree loaded from %s", filepath)
        return tree

    @staticmethod
    def majority_rule_consensus(trees: List[Tree], cutoff: float = 0.5) -> Tree:
        """Compute a majority-rule consensus tree.

        Args:
            trees: List of Bio.Phylo Tree objects.
            cutoff: Minimum clade frequency (0–1) to include in consensus.

        Returns:
            Consensus Tree object.
        """
        from Bio.Phylo.Consensus import majority_consensus

        consensus = majority_consensus(trees, cutoff=cutoff)
        logger.info("Computed majority-rule consensus from %d trees.", len(trees))
        return consensus

    @staticmethod
    def root_midpoint(tree: Tree) -> Tree:
        """Root the tree at the midpoint.

        Args:
            tree: Bio.Phylo Tree object.

        Returns:
            Rooted Tree.
        """
        tree.root_at_midpoint()
        return tree

    @staticmethod
    def terminal_names(tree: Tree) -> List[str]:
        """Return the names of all terminal (leaf) nodes.

        Args:
            tree: Bio.Phylo Tree object.

        Returns:
            List of terminal taxon names.
        """
        return [clade.name for clade in tree.get_terminals()]

    @staticmethod
    def star_tree(taxon_names: List[str]) -> Tree:
        """Create a star (fully unresolved) tree from a list of taxon names.

        Useful as a starting tree for heuristic searches.

        Args:
            taxon_names: List of taxon names.

        Returns:
            Bio.Phylo Tree with a single internal node connected to all taxa.
        """
        clades = [Clade(name=name) for name in taxon_names]
        root = Clade(clades=clades)
        return Tree(root=root, rooted=False)
