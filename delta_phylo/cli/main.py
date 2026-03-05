"""
CLI entry point for delta-phylo.

Run ``delta-phylo --help`` to see available commands.
"""

from __future__ import annotations

import click

from delta_phylo.cli.commands import (
    parse_cmd,
    matrix_cmd,
    parsimony_cmd,
    distance_cmd,
    export_cmd,
    plot_cmd,
    analyze_cmd,
)


@click.group()
@click.version_option(package_name="delta-phylo")
def cli() -> None:
    """delta-phylo: Morphological phylogenetic analyses using DELTA datasets.

    \b
    Workflow:
      1. delta-phylo parse    <dataset/>   - Inspect DELTA data
      2. delta-phylo matrix   <dataset/>   - Build morphological matrix
      3. delta-phylo parsimony <dataset/>  - Run parsimony analysis
      4. delta-phylo export nexus <dataset/> - Export to NEXUS
      5. delta-phylo plot tree.nwk         - Visualise tree
      6. delta-phylo analyze  <dataset/>   - Full pipeline (steps 1–5)
    """


cli.add_command(parse_cmd, name="parse")
cli.add_command(matrix_cmd, name="matrix")
cli.add_command(parsimony_cmd, name="parsimony")
cli.add_command(distance_cmd, name="distance")
cli.add_command(export_cmd, name="export")
cli.add_command(plot_cmd, name="plot")
cli.add_command(analyze_cmd, name="analyze")


if __name__ == "__main__":
    cli()
