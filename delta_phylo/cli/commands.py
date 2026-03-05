"""
CLI commands for delta-phylo.

Each command is a Click command group/command that calls the appropriate
delta-phylo module.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder
from delta_phylo.matrix.matrix_cleaning import MatrixCleaner
from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis
from delta_phylo.phylogeny.distance import DistanceAnalysis
from delta_phylo.phylogeny.tree_utils import TreeUtils
from delta_phylo.io.nexus_writer import NexusWriter
from delta_phylo.io.phylip_writer import PhylipWriter
from delta_phylo.io.tnt_writer import TNTWriter
from delta_phylo.io.csv_writer import CSVWriter
from delta_phylo.utils.validation import validate_delta_directory, validate_matrix
from delta_phylo.utils.logging import configure_root_logger

logger = logging.getLogger(__name__)


def _load_dataset(dataset_dir: str):
    """Helper: validate, read and return characters + taxa."""
    path = validate_delta_directory(dataset_dir)
    reader = DeltaReader(path)
    reader.read()
    return reader.characters, reader.taxa


def _build_matrix(dataset_dir: str):
    """Helper: load dataset and build MorphologicalMatrix."""
    chars, taxa = _load_dataset(dataset_dir)
    builder = MatrixBuilder(chars, taxa)
    matrix = builder.build()
    validate_matrix(matrix)
    return matrix


@click.command("parse")
@click.argument("dataset_dir")
def parse_cmd(dataset_dir: str) -> None:
    """Parse a DELTA dataset and display a summary.

    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    chars, taxa = _load_dataset(dataset_dir)
    click.echo(f"Characters : {len(chars)}")
    click.echo(f"Taxa       : {len(taxa)}")
    click.echo("\nFirst 5 characters:")
    for c in chars[:5]:
        click.echo(f"  {c}")
    click.echo("\nFirst 5 taxa:")
    for t in taxa[:5]:
        click.echo(f"  {t}")


@click.command("matrix")
@click.argument("dataset_dir")
@click.option("--output", "-o", default=None, help="Output CSV file path.")
@click.option(
    "--remove-constant",
    is_flag=True,
    default=False,
    help="Remove constant characters.",
)
def matrix_cmd(dataset_dir: str, output: str, remove_constant: bool) -> None:
    """Build and display the morphological matrix.

    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    matrix = _build_matrix(dataset_dir)
    if remove_constant:
        cleaner = MatrixCleaner(matrix)
        matrix = cleaner.remove_constant_characters()
    click.echo(matrix.df.to_string())
    if output:
        CSVWriter(matrix).write(output)
        click.echo(f"\nMatrix saved to {output}")


@click.command("parsimony")
@click.argument("dataset_dir")
@click.option("--replicates", "-r", default=10, show_default=True, type=int)
@click.option("--iterations", "-i", default=100, show_default=True, type=int)
@click.option("--output", "-o", default="parsimony_tree.nwk", show_default=True)
@click.option("--seed", default=42, show_default=True, type=int)
def parsimony_cmd(
    dataset_dir: str, replicates: int, iterations: int, output: str, seed: int
) -> None:
    """Run maximum parsimony analysis on a DELTA dataset.

    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    matrix = _build_matrix(dataset_dir)
    analysis = ParsimonyAnalysis(
        matrix=matrix,
        num_replicates=replicates,
        max_iterations=iterations,
        seed=seed,
    )
    trees = analysis.run()
    click.echo(f"Best parsimony score: {analysis.best_score}")
    click.echo(f"Number of equally parsimonious trees: {len(trees)}")
    best_tree = trees[0]
    TreeUtils.save_tree(best_tree, output)
    click.echo(f"Best tree saved to {output}")


@click.command("distance")
@click.argument("dataset_dir")
@click.option("--metric", "-m", default="hamming", show_default=True,
              type=click.Choice(["hamming", "simple_matching"]))
@click.option("--output", "-o", default="nj_tree.nwk", show_default=True)
def distance_cmd(dataset_dir: str, metric: str, output: str) -> None:
    """Run Neighbor Joining tree inference on a DELTA dataset.

    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    matrix = _build_matrix(dataset_dir)
    analysis = DistanceAnalysis(matrix=matrix, metric=metric)
    tree = analysis.run()
    TreeUtils.save_tree(tree, output)
    click.echo(f"NJ tree saved to {output}")


@click.command("export")
@click.argument("format", type=click.Choice(["nexus", "phylip", "tnt", "csv"]))
@click.argument("dataset_dir")
@click.option("--output", "-o", default=None, help="Output file path.")
def export_cmd(format: str, dataset_dir: str, output: str) -> None:
    """Export the morphological matrix to a phylogenetic file format.

    FORMAT: Output format (nexus, phylip, tnt, csv).
    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    matrix = _build_matrix(dataset_dir)
    extensions = {"nexus": ".nex", "phylip": ".phy", "tnt": ".tnt", "csv": ".csv"}
    if output is None:
        output = f"matrix{extensions[format]}"

    writers = {
        "nexus": NexusWriter,
        "phylip": PhylipWriter,
        "tnt": TNTWriter,
        "csv": CSVWriter,
    }
    writers[format](matrix).write(output)
    click.echo(f"Matrix exported to {output} ({format.upper()} format)")


@click.command("plot")
@click.argument("tree_file")
@click.option("--output", "-o", default=None, help="Output image file (PNG/PDF/SVG).")
@click.option("--title", default="Phylogenetic Tree", show_default=True)
@click.option("--format", "fmt", default="newick", show_default=True,
              type=click.Choice(["newick", "nexus"]))
def plot_cmd(tree_file: str, output: str, title: str, fmt: str) -> None:
    """Visualize a phylogenetic tree.

    TREE_FILE: Path to the tree file.
    """
    configure_root_logger()
    from delta_phylo.visualization.tree_plot import TreePlotter

    tree = TreeUtils.load_tree(tree_file, fmt=fmt)
    plotter = TreePlotter(tree)
    if output:
        plotter.save(output, title=title)
        click.echo(f"Tree plot saved to {output}")
    else:
        plotter.plot(title=title, show=False)
        click.echo("Tree plot rendered (no output file specified; use --output).")


@click.command("analyze")
@click.argument("dataset_dir")
@click.option("--method", "-m", default="parsimony",
              type=click.Choice(["parsimony", "nj"]), show_default=True)
@click.option("--output-dir", "-d", default=".", show_default=True)
@click.option("--replicates", "-r", default=10, show_default=True, type=int)
@click.option("--seed", default=42, show_default=True, type=int)
def analyze_cmd(
    dataset_dir: str,
    method: str,
    output_dir: str,
    replicates: int,
    seed: int,
) -> None:
    """Run the full analysis pipeline: parse → matrix → infer tree → export.

    DATASET_DIR: Path to the directory containing DELTA files.
    """
    configure_root_logger()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    click.echo("Step 1/4: Parsing DELTA dataset...")
    matrix = _build_matrix(dataset_dir)
    click.echo(f"  {matrix}")

    click.echo("Step 2/4: Exporting matrix...")
    NexusWriter(matrix).write(out / "matrix.nex")
    CSVWriter(matrix).write(out / "matrix.csv")

    click.echo(f"Step 3/4: Running {method} analysis...")
    if method == "parsimony":
        analysis = ParsimonyAnalysis(
            matrix=matrix, num_replicates=replicates, seed=seed
        )
        trees = analysis.run()
        best_tree = trees[0]
        tree_path = out / "best_tree.nwk"
        click.echo(f"  Best parsimony score: {analysis.best_score}")
    else:
        analysis_obj = DistanceAnalysis(matrix=matrix)
        best_tree = analysis_obj.run()
        tree_path = out / "nj_tree.nwk"

    TreeUtils.save_tree(best_tree, str(tree_path))
    click.echo(f"  Tree saved to {tree_path}")

    click.echo("Step 4/4: Generating tree plot...")
    from delta_phylo.visualization.tree_plot import TreePlotter

    plot_path = out / "tree.png"
    TreePlotter(best_tree).save(str(plot_path))
    click.echo(f"  Plot saved to {plot_path}")

    click.echo("\nAnalysis complete!")
    click.echo(f"Results in: {out.resolve()}")
