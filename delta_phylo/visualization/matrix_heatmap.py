"""
MatrixHeatmap: visualize morphological matrices as heatmaps using matplotlib.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


class MatrixHeatmap:
    """Visualise a morphological matrix as a heatmap.

    Args:
        matrix: MorphologicalMatrix to visualise.
    """

    def __init__(self, matrix) -> None:
        self.matrix = matrix

    def plot(
        self,
        title: str = "Morphological Matrix",
        figsize: Optional[Tuple[int, int]] = None,
        output_file: Optional[str | Path] = None,
        show: bool = True,
        cmap: str = "viridis",
    ) -> None:
        """Render the morphological matrix as a heatmap.

        Missing values are shown in grey.

        Args:
            title: Figure title.
            figsize: Figure dimensions (width, height) in inches.
                If None, auto-scales based on matrix size.
            output_file: Path to save the figure. Format inferred from extension.
            show: Whether to call ``plt.show()``.
            cmap: Matplotlib colormap name.
        """
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError as exc:
            raise ImportError(
                "matplotlib is required for heatmap plots. "
                "Install with: pip install matplotlib"
            ) from exc

        df = self.matrix.df
        n_taxa, n_chars = df.shape

        if figsize is None:
            w = max(6, min(n_chars // 3, 40))
            h = max(4, min(n_taxa // 2, 30))
            figsize = (w, h)

        fig, ax = plt.subplots(figsize=figsize)

        data = df.values.astype(float)
        # NaN → masked array for proper colouring
        masked = np.ma.masked_invalid(data)
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad(color="lightgrey")

        im = ax.imshow(masked, aspect="auto", cmap=cmap_obj, interpolation="nearest")
        ax.set_title(title)
        ax.set_xlabel("Characters")
        ax.set_ylabel("Taxa")

        # Axis labels
        if n_chars <= 50:
            ax.set_xticks(range(n_chars))
            ax.set_xticklabels(df.columns, rotation=90, fontsize=7)
        if n_taxa <= 50:
            ax.set_yticks(range(n_taxa))
            ax.set_yticklabels(df.index, fontsize=7)

        plt.colorbar(im, ax=ax, label="State")

        if output_file is not None:
            output_file = Path(output_file)
            fmt = output_file.suffix.lstrip(".").lower() or "png"
            fig.savefig(output_file, format=fmt, bbox_inches="tight", dpi=150)
            logger.info("Matrix heatmap saved to %s", output_file)

        if show:
            plt.show()

        plt.close(fig)

    def save(
        self,
        output_file: str | Path,
        title: str = "Morphological Matrix",
        figsize: Optional[Tuple[int, int]] = None,
    ) -> None:
        """Save the heatmap without displaying it.

        Args:
            output_file: Output file path.
            title: Figure title.
            figsize: Figure dimensions in inches.
        """
        self.plot(title=title, figsize=figsize, output_file=output_file, show=False)
