# delta-phylo

**delta-phylo** is a modern, open-source Python library for morphological
phylogenetic analyses using DELTA format datasets. It is designed for
systematists, taxonomists, and evolutionary biologists working with
morphological characters.

---

## Features

- **DELTA Parser** – Load and parse standard DELTA files (`characters`, `items`, `specs`)
  - Ordered multistate characters (`OM`/`RU` directive → `CharacterType.ORDERED`)
  - Numeric/continuous characters (`RN` directive → `CharacterType.NUMERIC`)
  - Polymorphic scores, ranges, and missing/inapplicable data
- **Morphological Matrix** – Convert DELTA data to structured NumPy/Pandas matrices
  - Polymorphic resolution strategy: `"first"` (min state), `"missing"` (NaN), or `"majority"` (mode)
  - `to_numpy_float()` preserves decimal values for numeric characters
- **Matrix Cleaning**
  - `remove_constant_characters()` – drops invariant columns
  - `remove_autapomorphies()` – removes characters where every derived state is unique to a single taxon; synapomorphies shared by ≥ 2 taxa are preserved
  - `filter_by_completeness()` – filters taxa/characters below a completeness threshold
  - `impute_missing()` – mode or zero imputation
- **Multiple phylogenetic methods**:
  - Maximum Parsimony (Fitch algorithm + NNI heuristic search)
  - Distance-based Neighbor Joining – Hamming, Simple Matching, and **Gower distance** (handles mixed discrete + numeric data automatically)
  - Maximum Likelihood – **Felsenstein pruning** under the Mk model (topology- and branch-length-aware)
- **Matrix Export**:
  - NEXUS (sequential and **interleaved** format with configurable `block_width`)
  - PHYLIP, TNT, CSV
- **Tree Visualization** – Static plots (PNG, PDF, SVG) via matplotlib
- **Matrix Heatmap** – Visualize completeness and character distribution
- **CLI** – Command-line interface for all operations
- **Reproducible** – Seed-controlled random analyses
- **Full public API** – Every major class and function importable directly from `delta_phylo`

---

## Installation

### From source (development)

```bash
git clone https://github.com/nathanufpb/morfology_analysis.git
cd morfology_analysis
pip install -e ".[dev]"
```

### From PyPI (when published)

```bash
pip install delta-phylo
```

---

## Quickstart

### Python API

```python
from delta_phylo import (
    DeltaReader, MatrixBuilder, MatrixCleaner,
    DistanceAnalysis, ParsimonyAnalysis, LikelihoodAnalysis,
    NexusWriter, TreeUtils,
)

# 1. Parse DELTA dataset
reader = DeltaReader("examples/example_delta_dataset")
reader.read()

# 2. Build matrix (polymorphic taxa resolved by majority-vote)
matrix = MatrixBuilder(reader.characters, reader.taxa,
                       polymorphic_strategy="majority").build()
print(matrix.df)

# 3. Clean matrix
cleaner = MatrixCleaner(matrix)
matrix = cleaner.remove_constant_characters()
matrix = cleaner.remove_autapomorphies()

# 4. Infer Neighbor Joining tree (Gower distance used automatically
#    when the matrix contains numeric characters)
nj = DistanceAnalysis(matrix=matrix, metric="hamming")
tree = nj.run()
print(TreeUtils.tree_to_newick(tree))

# 5. Score tree with Maximum Likelihood (Felsenstein pruning, Mk model)
ll = LikelihoodAnalysis(matrix).score_tree(tree)
print(f"Log-likelihood: {ll:.4f}")

# 6. Export to NEXUS – interleaved blocks of 60 characters
NexusWriter(matrix).write("matrix.nex", interleave=True, block_width=60)
```

### Command Line

```bash
# Parse dataset
delta-phylo parse examples/example_delta_dataset/

# Build and display matrix
delta-phylo matrix examples/example_delta_dataset/

# Run parsimony analysis
delta-phylo parsimony examples/example_delta_dataset/ --output mp_tree.nwk

# Run Neighbor Joining
delta-phylo distance examples/example_delta_dataset/ --output nj_tree.nwk

# Export to NEXUS
delta-phylo export nexus examples/example_delta_dataset/

# Visualize tree
delta-phylo plot mp_tree.nwk --output tree.png

# Full pipeline
delta-phylo analyze examples/example_delta_dataset/ --output-dir results/
```

---

## Numeric / Continuous Characters

Characters declared with `*CHARACTER TYPES RN <ids>` in the DELTA `specs` file
are parsed as `CharacterType.NUMERIC`. Their raw float values are stored in the
matrix and preserved by `MorphologicalMatrix.to_numpy_float()`.

When a matrix contains numeric characters, `DistanceAnalysis` automatically
switches to **Gower's distance** — a mixed-data metric that normalises numeric
differences by per-column value range and uses 0/1 mismatch for discrete
characters. You can also request it explicitly:

```python
from delta_phylo import DistanceAnalysis, gower_distance

analysis = DistanceAnalysis(matrix, metric="gower")
tree = analysis.run()
```

---

## NEXUS Interleaved Format

```python
from delta_phylo import NexusWriter

# Sequential (default)
NexusWriter(matrix).write("matrix.nex")

# Interleaved – adds INTERLEAVE=YES to FORMAT, splits MATRIX into blocks
NexusWriter(matrix).write("matrix.nex", interleave=True, block_width=60)
```

---

## Project Architecture

```
delta-phylo
│
├── delta_phylo/
│   ├── parser/          DELTA file parsing
│   │   ├── delta_reader.py   DeltaReader (ordered + numeric char detection)
│   │   ├── characters.py     Character & CharacterType
│   │   │                     (UNORDERED, ORDERED, BINARY, NUMERIC)
│   │   ├── taxa.py           Taxon class
│   │   └── states.py         CharacterState class
│   │
│   ├── matrix/          Morphological matrix construction
│   │   ├── matrix_builder.py   MatrixBuilder → MorphologicalMatrix
│   │   │                       (strategies: first / missing / majority)
│   │   ├── matrix_cleaning.py  MatrixCleaner (autapomorphies, completeness)
│   │   └── encoding.py         MatrixEncoder (binary encoding, symbol strings)
│   │
│   ├── phylogeny/       Phylogenetic analysis
│   │   ├── parsimony.py    Fitch + NNI parsimony search
│   │   ├── distance.py     Hamming / SMC / Gower + Neighbor Joining
│   │   ├── likelihood.py   Felsenstein pruning, Mk model
│   │   └── tree_utils.py   Tree I/O, consensus, utilities
│   │
│   ├── io/              Matrix/tree export
│   │   ├── nexus_writer.py   Sequential & interleaved NEXUS
│   │   ├── phylip_writer.py
│   │   ├── tnt_writer.py
│   │   └── csv_writer.py
│   │
│   ├── visualization/   Plotting
│   │   ├── tree_plot.py       TreePlotter
│   │   └── matrix_heatmap.py  MatrixHeatmap
│   │
│   ├── cli/             Command line interface
│   │   ├── main.py      Entry point (delta-phylo)
│   │   └── commands.py  Individual commands
│   │
│   └── utils/           Utilities
│       ├── validation.py
│       └── logging.py
│
├── tests/               Automated pytest tests (220 tests, 0 failures)
├── examples/            Example DELTA dataset + run_analysis.py
└── docs/                Documentation
```

---

## Running the Example

```bash
python examples/run_analysis.py
```

Output files will be written to `examples/output/`.

---

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

Current status: **220 tests, 0 failures**.

---

## Dependencies

| Package    | Purpose                                      |
|------------|----------------------------------------------|
| numpy      | Efficient array operations                   |
| pandas     | DataFrame-based matrix storage               |
| biopython  | Tree structures and I/O                      |
| matplotlib | Tree and matrix visualization                |
| scipy      | Distance computations                        |
| click      | Command line interface                       |

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citation

If you use delta-phylo in your research, please cite:

> delta-phylo contributors (2024). *delta-phylo: Morphological phylogenetic
> analyses using DELTA format datasets*. https://github.com/nathanufpb/morfology_analysis
