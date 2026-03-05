"""
TreePlotter: static tree visualization using matplotlib and Bio.Phylo.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


class TreePlotter:
    """Plot phylogenetic trees using matplotlib.

    Args:
        tree: Bio.Phylo Tree object to visualise.
    """

    def __init__(self, tree) -> None:
        self.tree = tree

    def plot(
        self,
        title: str = "Phylogenetic Tree",
        figsize: Tuple[int, int] = (10, 8),
        output_file: Optional[str | Path] = None,
        show: bool = True,
    ) -> None:
        """Render the tree as a matplotlib figure.

        Args:
            title: Figure title.
            figsize: Matplotlib figure size (width, height) in inches.
            output_file: If provided, save the figure to this path.
                Supported formats: PNG, PDF, SVG (inferred from extension).
            show: Whether to call ``plt.show()`` after rendering.
        """
        try:
            import matplotlib
            matplotlib.use("Agg")  # Non-interactive backend
            import matplotlib.pyplot as plt
            from Bio import Phylo
        except ImportError as exc:
            raise ImportError(
                "matplotlib and biopython are required for tree plotting. "
                "Install with: pip install matplotlib biopython"
            ) from exc

        fig, ax = plt.subplots(figsize=figsize)
        Phylo.draw(self.tree, axes=ax, do_show=False)
        ax.set_title(title)

        if output_file is not None:
            output_file = Path(output_file)
            fmt = output_file.suffix.lstrip(".").lower() or "png"
            fig.savefig(output_file, format=fmt, bbox_inches="tight", dpi=150)
            logger.info("Tree plot saved to %s", output_file)

        if show:
            plt.show()

        plt.close(fig)

    def save(
        self,
        output_file: str | Path,
        title: str = "Phylogenetic Tree",
        figsize: Tuple[int, int] = (10, 8),
    ) -> None:
        """Save the tree figure without displaying it.

        Args:
            output_file: Output file path (PNG, PDF, or SVG).
            title: Figure title.
            figsize: Figure dimensions in inches.
        """
        self.plot(title=title, figsize=figsize, output_file=output_file, show=False)
