"""
Example script demonstrating the full delta-phylo analysis pipeline.

Usage::

    python examples/run_analysis.py

This script performs:
  1. Parse DELTA dataset
  2. Build morphological matrix
  3. Export matrix to NEXUS, PHYLIP, TNT, and CSV
  4. Run Neighbor Joining tree inference
  5. Run Maximum Parsimony tree inference
  6. Generate tree and matrix visualizations
"""

from __future__ import annotations

import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run the complete analysis pipeline."""
    # ---------------------------------------------------------------------------
    # Paths
    # ---------------------------------------------------------------------------
    examples_dir = Path(__file__).parent
    dataset_dir = examples_dir / "example_delta_dataset"
    output_dir = examples_dir / "output"
    output_dir.mkdir(exist_ok=True)

    # ---------------------------------------------------------------------------
    # Step 1: Parse DELTA dataset
    # ---------------------------------------------------------------------------
    logger.info("=== Step 1: Parsing DELTA dataset ===")
    from delta_phylo.parser.delta_reader import DeltaReader

    reader = DeltaReader(dataset_dir)
    reader.read()
    logger.info("Parsed %d characters and %d taxa.", len(reader.characters), len(reader.taxa))

    # ---------------------------------------------------------------------------
    # Step 2: Build morphological matrix
    # ---------------------------------------------------------------------------
    logger.info("=== Step 2: Building morphological matrix ===")
    from delta_phylo.matrix.matrix_builder import MatrixBuilder
    from delta_phylo.matrix.matrix_cleaning import MatrixCleaner

    builder = MatrixBuilder(reader.characters, reader.taxa)
    matrix = builder.build()
    logger.info("Matrix shape: %s", matrix.shape)
    print("\nMorphological matrix (first 5 columns):")
    print(matrix.df.iloc[:, :5])

    # ---------------------------------------------------------------------------
    # Step 3: Export matrix
    # ---------------------------------------------------------------------------
    logger.info("=== Step 3: Exporting matrix ===")
    from delta_phylo.io.nexus_writer import NexusWriter
    from delta_phylo.io.phylip_writer import PhylipWriter
    from delta_phylo.io.tnt_writer import TNTWriter
    from delta_phylo.io.csv_writer import CSVWriter

    NexusWriter(matrix).write(output_dir / "matrix.nex")
    PhylipWriter(matrix).write(output_dir / "matrix.phy")
    TNTWriter(matrix).write(output_dir / "matrix.tnt")
    CSVWriter(matrix).write(output_dir / "matrix.csv")
    logger.info("Matrix exported to %s", output_dir)

    # ---------------------------------------------------------------------------
    # Step 4: Neighbor Joining tree
    # ---------------------------------------------------------------------------
    logger.info("=== Step 4: Neighbor Joining analysis ===")
    from delta_phylo.phylogeny.distance import DistanceAnalysis
    from delta_phylo.phylogeny.tree_utils import TreeUtils

    nj_analysis = DistanceAnalysis(matrix=matrix, metric="hamming")
    nj_tree = nj_analysis.run()
    nj_path = output_dir / "nj_tree.nwk"
    TreeUtils.save_tree(nj_tree, str(nj_path))
    logger.info("NJ tree saved to %s", nj_path)
    print("\nNJ tree (Newick):")
    print(TreeUtils.tree_to_newick(nj_tree))

    # ---------------------------------------------------------------------------
    # Step 5: Maximum Parsimony analysis
    # ---------------------------------------------------------------------------
    logger.info("=== Step 5: Maximum Parsimony analysis ===")
    from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis

    parsimony = ParsimonyAnalysis(matrix=matrix, num_replicates=5, seed=42)
    mp_trees = parsimony.run()
    logger.info("Best parsimony score: %d", parsimony.best_score)
    logger.info("Equally parsimonious trees: %d", len(mp_trees))
    mp_path = output_dir / "parsimony_tree.nwk"
    TreeUtils.save_tree(mp_trees[0], str(mp_path))
    logger.info("Best MP tree saved to %s", mp_path)

    # ---------------------------------------------------------------------------
    # Step 6: Visualizations
    # ---------------------------------------------------------------------------
    logger.info("=== Step 6: Generating visualizations ===")
    from delta_phylo.visualization.tree_plot import TreePlotter
    from delta_phylo.visualization.matrix_heatmap import MatrixHeatmap

    TreePlotter(nj_tree).save(
        str(output_dir / "nj_tree.png"), title="Neighbor Joining Tree"
    )
    TreePlotter(mp_trees[0]).save(
        str(output_dir / "parsimony_tree.png"), title="Maximum Parsimony Tree"
    )
    MatrixHeatmap(matrix).save(
        str(output_dir / "matrix_heatmap.png"), title="Morphological Matrix"
    )

    logger.info("=== Analysis complete! ===")
    logger.info("All results saved in: %s", output_dir.resolve())
    print(f"\nResults saved in: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
