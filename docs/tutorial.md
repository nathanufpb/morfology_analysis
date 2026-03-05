# Tutorial: Morphological Phylogenetics with delta-phylo

This tutorial walks through a complete morphological phylogenetic analysis
using the example arthropod dataset included with delta-phylo.

---

## Step 0: Install delta-phylo

```bash
git clone https://github.com/nathanufpb/morfology_analysis.git
cd morfology_analysis
pip install -e ".[dev]"
```

Verify the installation:

```bash
delta-phylo --version
delta-phylo --help
```

---

## Step 1: Inspect the DELTA dataset

The example dataset contains 10 arthropod taxa and 20 morphological characters.

```bash
delta-phylo parse examples/example_delta_dataset/
```

Expected output:

```
Characters : 20
Taxa       : 10

First 5 characters:
  Character(id=1, name='body size', type=UNORDERED, states=3)
  Character(id=2, name='body symmetry', type=BINARY, states=2)
  ...
```

---

## Step 2: Build the morphological matrix

```bash
delta-phylo matrix examples/example_delta_dataset/ --output matrix.csv
```

Or in Python:

```python
from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.matrix.matrix_builder import MatrixBuilder

reader = DeltaReader("examples/example_delta_dataset")
reader.read()
builder = MatrixBuilder(reader.characters, reader.taxa)
matrix = builder.build()
print(matrix.df)
```

The matrix is a `taxa × characters` table where each cell is a 0-based
integer state, or `NaN` for missing data.

---

## Step 3: Clean the matrix

Remove uninformative characters:

```python
from delta_phylo.matrix.matrix_cleaning import MatrixCleaner

cleaner = MatrixCleaner(matrix)
matrix = cleaner.remove_constant_characters()
print(f"After cleaning: {matrix.shape}")
```

---

## Step 4: Export the matrix

```bash
# NEXUS (for MrBayes, PAUP*)
delta-phylo export nexus examples/example_delta_dataset/ --output matrix.nex

# PHYLIP (for PHYLIP package)
delta-phylo export phylip examples/example_delta_dataset/ --output matrix.phy

# TNT (for TNT parsimony)
delta-phylo export tnt examples/example_delta_dataset/ --output matrix.tnt
```

---

## Step 5: Infer a Neighbor Joining tree

```bash
delta-phylo distance examples/example_delta_dataset/ --output nj_tree.nwk
```

Or in Python:

```python
from delta_phylo.phylogeny.distance import DistanceAnalysis
from delta_phylo.phylogeny.tree_utils import TreeUtils

analysis = DistanceAnalysis(matrix=matrix, metric="hamming")
tree = analysis.run()
print(TreeUtils.tree_to_newick(tree))
TreeUtils.save_tree(tree, "nj_tree.nwk")
```

---

## Step 6: Run Maximum Parsimony

```bash
delta-phylo parsimony examples/example_delta_dataset/ \
    --replicates 20 --iterations 200 --output mp_tree.nwk
```

Or in Python:

```python
from delta_phylo.phylogeny.parsimony import ParsimonyAnalysis

analysis = ParsimonyAnalysis(matrix=matrix, num_replicates=20, seed=42)
trees = analysis.run()
print(f"Best score: {analysis.best_score}")
print(f"Trees found: {len(trees)}")

# Consensus tree
consensus = analysis.consensus_tree(cutoff=0.5)
TreeUtils.save_tree(consensus, "consensus_tree.nwk")
```

---

## Step 7: Visualize the tree

```bash
delta-phylo plot mp_tree.nwk --output tree.png --title "Arthropod Phylogeny"
```

Or in Python:

```python
from delta_phylo.visualization.tree_plot import TreePlotter

plotter = TreePlotter(tree)
plotter.save("tree.png", title="Arthropod Phylogeny (NJ)")
```

---

## Step 8: Visualize the matrix as a heatmap

```python
from delta_phylo.visualization.matrix_heatmap import MatrixHeatmap

heatmap = MatrixHeatmap(matrix)
heatmap.save("matrix_heatmap.png", title="Morphological Data Completeness")
```

---

## Full Pipeline (one command)

```bash
delta-phylo analyze examples/example_delta_dataset/ \
    --method parsimony \
    --output-dir results/ \
    --replicates 10
```

This runs the complete workflow:

1. Parse DELTA dataset
2. Build matrix
3. Export NEXUS and CSV
4. Run parsimony analysis
5. Save tree and visualization

---

## Tips

- For large datasets (>100 taxa), use `--replicates 50` or more.
- Export to NEXUS and run **MrBayes** for Bayesian inference with the Mk model.
- Export to NEXUS and run **IQ-TREE** with `-m MK+G` for ML inference.
- Use `--remove-constant` with the `matrix` command to remove invariant characters.
- Use `MatrixCleaner.filter_by_completeness()` to remove highly incomplete data.
