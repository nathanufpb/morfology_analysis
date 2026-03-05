# Architecture

## Overview

delta-phylo follows a modular, layered architecture designed for scientific
reproducibility and extensibility.

```
DELTA Dataset (files)
        │
        ▼
  ┌─────────────┐
  │   Parser    │  delta_phylo/parser/
  │  (DeltaReader, Character, Taxon, CharacterState)
  └──────┬──────┘
         │ List[Character], List[Taxon]
         ▼
  ┌─────────────┐
  │   Matrix    │  delta_phylo/matrix/
  │  (MatrixBuilder → MorphologicalMatrix)
  └──────┬──────┘
         │ MorphologicalMatrix (pandas DataFrame + metadata)
         ├──────────────────────────────────┐
         ▼                                  ▼
  ┌────────────┐                    ┌───────────────┐
  │  Phylogeny │                    │      IO       │
  │ (Parsimony,│                    │ (NEXUS, PHYLIP│
  │  Distance, │                    │  TNT, CSV)    │
  │  TreeUtils)│                    └───────────────┘
  └──────┬─────┘
         │ Bio.Phylo Tree
         ▼
  ┌──────────────┐
  │Visualization │  delta_phylo/visualization/
  │(TreePlotter, │
  │ MatrixHeatmap│
  └──────────────┘
```

## Data Flow

1. **Parse** → `DeltaReader` reads DELTA files into `Character`, `Taxon`, and
   `CharacterState` objects.
2. **Build** → `MatrixBuilder` converts the parsed objects into a
   `MorphologicalMatrix` (a pandas DataFrame with metadata).
3. **Clean** → `MatrixCleaner` removes uninformative characters and filters
   by completeness.
4. **Encode** → `MatrixEncoder` converts multistate characters to binary or
   symbol strings for export.
5. **Analyse** → `ParsimonyAnalysis` or `DistanceAnalysis` infers a
   phylogenetic tree (Bio.Phylo Tree object).
6. **Export** → IO writers serialize the matrix to NEXUS, PHYLIP, TNT, or CSV.
7. **Visualize** → `TreePlotter` and `MatrixHeatmap` create matplotlib figures.

## Key Classes

| Class | Module | Responsibility |
|-------|--------|---------------|
| `DeltaReader` | `parser.delta_reader` | Parse DELTA files |
| `Character` | `parser.characters` | Store character definition |
| `Taxon` | `parser.taxa` | Store taxon + scores |
| `CharacterState` | `parser.states` | Individual state |
| `MorphologicalMatrix` | `matrix.matrix_builder` | Taxa × characters DataFrame |
| `MatrixBuilder` | `matrix.matrix_builder` | Build matrix from parsed data |
| `MatrixCleaner` | `matrix.matrix_cleaning` | Filter/clean matrix |
| `MatrixEncoder` | `matrix.encoding` | Encode for export |
| `ParsimonyAnalysis` | `phylogeny.parsimony` | MP inference |
| `DistanceAnalysis` | `phylogeny.distance` | NJ inference |
| `TreeUtils` | `phylogeny.tree_utils` | Tree manipulation |
| `NexusWriter` | `io.nexus_writer` | NEXUS export |
| `PhylipWriter` | `io.phylip_writer` | PHYLIP export |
| `TNTWriter` | `io.tnt_writer` | TNT export |
| `CSVWriter` | `io.csv_writer` | CSV export |
| `TreePlotter` | `visualization.tree_plot` | Tree plotting |
| `MatrixHeatmap` | `visualization.matrix_heatmap` | Heatmap plotting |

## Extensibility

To add a new phylogenetic method:
1. Create a new module in `delta_phylo/phylogeny/`.
2. Accept a `MorphologicalMatrix` and return a `Bio.Phylo Tree`.
3. Register a CLI command in `delta_phylo/cli/commands.py`.

To add a new export format:
1. Create a new writer in `delta_phylo/io/`.
2. Implement a `write(filepath)` method.
3. Register in the CLI `export` command.
