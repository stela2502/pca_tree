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

---

## Installation

`pca_tree` depends on openblas and we therefore need to be able to compile this

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
## Using from R

You can call the pca_tree Rust binary directly from R using this helper function:
```R
pca_tree <- function(
    seqs,
    k = 3,
    prefix = tempfile("pca_tree_"),
    bin = "pca_tree",
    plots = FALSE
) {
    # Ensure output directory exists
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

    # Derive file paths
    infile   <- paste0(prefix, "_input.txt")
    coords   <- paste0(prefix, "_pca.tsv")
    edges    <- paste0(prefix, "_tree.tsv")
    pca_png  <- paste0(prefix, "_pca.png")
    tree_png <- paste0(prefix, "_tree.png")

    # Write sequences
    writeLines(seqs, infile)

    # Build CLI args
    args <- c(
        infile,
        "--k", k,
        "--coords", coords,
        "--edges", edges
    )

    if (plots) {
        args <- c(args,
                  "--plot-pca", pca_png,
                  "--plot-tree", tree_png)
    }

    # Run Rust binary
    out <- system2(bin, args, stdout = TRUE, stderr = TRUE)

    # Load results
    coords_mat <- as.matrix(read.table(coords))
    edges_df   <- read.table(edges,
                             col.names = c("parent", "child", "dist"))

    list(
        coords = coords_mat,
        edges  = edges_df,
        prefix = prefix,
        pca_plot  = if (plots) pca_png else NULL,
        tree_plot = if (plots) tree_png else NULL
    )
}
```
So assuming you have been using Change-O to analyze VDJ recombination evens you wil end up at something like "YourSample_ProductiveCloneDefined.tsv" or something like that.
These files would contain an sequence_alignment column and a column_id column and with that you can uste the above function like this:

```R
my_infile="Path/To/Your/file.tsv"
df = read_delim( my_infile )
clones = tables(df[,'clone_id'])
## take the Change-O for the smaller if you want
clones = names(clones[which(clones > 200)])
for (clone_id in clones){
    srts = df[which(df[,'clone_id'] == clone_id), "sequence_alignment"]
    pca_tree( srts, k=30, prefix= basename(my_infile), plots=T)
}

```

And now you have a lot of PCS coordinates, PCA  and Tree plots.
Not that I have done anything with them at this time


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
