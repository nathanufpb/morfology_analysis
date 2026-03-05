# delta-phylo

**delta-phylo** is a modern, open-source Python library for morphological
phylogenetic analyses using DELTA format datasets. It is designed for
systematists, taxonomists, and evolutionary biologists working with
morphological characters.

---

## Features

- **DELTA Parser** – Load and parse standard DELTA files (`characters`, `items`, `specs`)
- **Morphological Matrix** – Convert DELTA data into structured NumPy/Pandas matrices
- **Multiple phylogenetic methods**:
  - Maximum Parsimony (Fitch algorithm + NNI heuristic search)
  - Distance-based Neighbor Joining (Hamming & Simple Matching)
  - Maximum Likelihood (simplified Mk model)
- **Matrix Export** – NEXUS, PHYLIP, TNT, CSV formats
- **Tree Visualization** – Static plots (PNG, PDF, SVG) using matplotlib
- **Matrix Heatmap** – Visualize completeness and character distribution
- **CLI** – Simple command-line interface for all operations
- **Reproducible** – Seed-controlled random analyses

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
from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder
from delta_phylo.phylogeny.distance import DistanceAnalysis
from delta_phylo.phylogeny.tree_utils import TreeUtils
from delta_phylo.io.nexus_writer import NexusWriter

# 1. Parse DELTA dataset
reader = DeltaReader("examples/example_delta_dataset")
reader.read()

# 2. Build matrix
builder = MatrixBuilder(reader.characters, reader.taxa)
matrix = builder.build()
print(matrix.df)

# 3. Infer Neighbor Joining tree
nj = DistanceAnalysis(matrix=matrix, metric="hamming")
tree = nj.run()
print(TreeUtils.tree_to_newick(tree))

# 4. Export to NEXUS
NexusWriter(matrix).write("matrix.nex")
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

## Project Architecture

```
delta-phylo
│
├── delta_phylo/
│   ├── parser/          DELTA file parsing
│   │   ├── delta_reader.py   Main reader (DeltaReader)
│   │   ├── characters.py     Character & CharacterType classes
│   │   ├── taxa.py           Taxon class
│   │   └── states.py         CharacterState class
│   │
│   ├── matrix/          Morphological matrix construction
│   │   ├── matrix_builder.py   MatrixBuilder → MorphologicalMatrix
│   │   ├── matrix_cleaning.py  MatrixCleaner (constant chars, completeness)
│   │   └── encoding.py         MatrixEncoder (binary, symbol strings)
│   │
│   ├── phylogeny/       Phylogenetic analysis
│   │   ├── parsimony.py    Fitch + NNI parsimony search
│   │   ├── distance.py     Hamming/SMC + Neighbor Joining
│   │   ├── likelihood.py   Simplified Mk model
│   │   └── tree_utils.py   Tree I/O, consensus, utilities
│   │
│   ├── io/              Matrix/tree export
│   │   ├── nexus_writer.py
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
├── tests/               Automated pytest tests
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

---

## Dependencies

| Package    | Purpose                        |
|------------|--------------------------------|
| numpy      | Efficient array operations     |
| pandas     | DataFrame-based matrix storage |
| biopython  | Tree structures and I/O        |
| matplotlib | Tree and matrix visualization  |
| scipy      | Distance computations          |
| click      | Command line interface         |

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Citation

If you use delta-phylo in your research, please cite:

> delta-phylo contributors (2024). *delta-phylo: Morphological phylogenetic
> analyses using DELTA format datasets*. https://github.com/nathanufpb/morfology_analysis