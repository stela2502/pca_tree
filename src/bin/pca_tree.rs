use clap::Parser;
use pca_tree::PcaTree;
use std::fs::read_to_string;
use std::path::{Path, PathBuf};


#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    /// Path to file with one DNA sequence per line
    input: String,

    /// Number of PCA components
    #[arg(short, long, default_value_t = 3)]
    k: usize,

    /// Output table with PCA coordinates
    #[arg(long)]
    coords: Option<String>,

    /// Output table with MST edges
    #[arg(long)]
    edges: Option<String>,

    /// Write PCA plot (PNG)
    #[arg(long)]
    plot_pca: Option<String>,

    /// Write PCA tree plot (PNG)
    #[arg(long)]
    plot_tree: Option<String>,

}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Args::parse();   // ✅ parse ONCE

    let input = read_to_string(&args.input)?;

    let seqs: Vec<String> = input
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|s| s.to_string())
        .collect();

    // ✅ NEW: constructor does all work
    let model = PcaTree::new(seqs, args.k)?;


    // Determine PCA output path
    let coords_path: PathBuf = if let Some(user) = args.coords.as_ref() {
        PathBuf::from(user)
    } else {
        default_output_path(&args.input, "_pca.tsv")
    };

    // Determine MST edge output path
    let edges_path: PathBuf = if let Some(user) = args.edges.as_ref() {
        PathBuf::from(user)
    } else {
        default_output_path(&args.input, "_tree.tsv")
    };

    model.pca.to_tsv(&coords_path)?;
    println!("Written PCA coords → {}", coords_path.display());
    
    model.tree.to_tsv(&edges_path)?;
    println!("Written MSt edges → {}", edges_path.display());
    

    #[cfg(feature = "plot")]
    {
        if let Some(f) = args.plot_pca {
            model.pca.plot_2d_clusters( &model.tree, &f)?;
            eprintln!("✅ PCA plot written to {f}");
        }

        if let Some(f) = args.plot_tree {
            model.tree.plot_2d(&model.coords(), &f)?;
            eprintln!("✅ Tree plot written to {f}");
        }
    }


    #[cfg(not(feature = "plot"))]
    {
        if args.plot_pca.is_some() || args.plot_tree.is_some() {
            eprintln!("⚠️ Plotting is disabled. Recompile with: cargo build --features plot");
        }
    }
    

    Ok(())
}



fn default_output_path(input: &str, suffix: &str) -> PathBuf {
    let path = Path::new(input);
    let stem = path.file_stem().unwrap_or_default();
    let parent = path.parent().unwrap_or_else(|| Path::new("."));

    let mut out = PathBuf::new();
    out.push(parent);
    out.push(format!("{}{}", stem.to_string_lossy(), suffix));

    out
}

