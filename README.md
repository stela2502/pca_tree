# pca_tree

High-performance PCA + MST engine for large-scale BCR / immunoglobulin sequence analysis.

`pca_tree` is a Rust library and CLI tool designed for extremely fast geometric
exploration of antibody repertoires. It provides:

- sequence encoding (one-hot or consensus-relative)
- PCA dimensionality reduction
- Minimum-Spanning-Tree (MST) construction
- clustering heuristics
- optional plotting using `plotters`
- bindings for R via extendr (see pcaTreeR)

## Features

### 🚀 Fast PCA
Implemented using `ndarray` and optional BLAS acceleration.

### 🌳 MST Lineage Geometry
Computes a geometric minimum-spanning-tree over PCA embeddings, enabling:
- clone structure detection
- lineage shape exploration
- cluster segmentation
- downstream integration with IgPhyML or Change-O

### 🎨 Optional Plotting
Produces PNG visualizations of:
- PCA scatter plots
- MST-overlaid PCA plots
using the `plotters` crate.

### 🔗 R Integration
An accompanying R package (`pcaTreeR`) provides high-level access to this library.

---

## Installation

### From source:

```bash
git clone https://github.com/stela2502/pca_tree.git
cd pca_tree
cargo build --release
```

### Install CLI globally:

```bash
cargo install --path .
```

---

## CLI Usage

```bash
pca_tree <sequences.txt> -k 3
```

Input format:
- One DNA sequence per line
- All sequences must be same length
- Ambiguous bases resolved by consensus encoding

Example output:

```
PCA coordinates:
0   [0.12, -0.04, 0.88]
1   [0.11, -0.05, 0.91]
...

Tree edges (parent child dist):
9337    1   0.0012
2508    4   0.0711
...
```

---

## Library Usage

Add to `Cargo.toml`:

```toml
[dependencies]
pca_tree = { git = "https://github.com/stela2502/pca_tree.git" }
```

Example:

```rust
use pca_tree::PcaTree;

let seqs = vec!["ACGTACGT".into(), "ACGTTCGT".into()];
let mut model = PcaTree::new(seqs, 3);
model.fit().unwrap();

println!("{:?}", model.coords());
println!("{:?}", model.tree());
```

---

## Plotting

```rust
model.pca.plot_2d("pca.png")?;
model.tree.plot_2d(model.coords,  "mst.png")?;
```

---

## Performance Notes

- Easily handles **10,000+ sequences**
- Pure Rust MST and PCA steps
- No recursion or deeply nested data structures
- Scales linearly in memory and time
- Compatible with musl, HPC, and bindgen

---

## License

MIT License
