use clap::Parser;
use pca_tree::PcaTree;
use std::fs::read_to_string;

#[derive(Parser)]
#[command(author, version, about)]
struct Args {
    /// Path to file with one DNA sequence per line
    input: String,

    /// Number of PCA components
    #[arg(short, long, default_value_t = 3)]
    k: usize,

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

    println!("PCA coordinates:");
    for (i, row) in model.coords().outer_iter().enumerate() {
        println!("{i}\t{:?}", row.to_vec());
    }

    println!("\nTree edges (parent child dist):");
    for (p, c, d) in model.tree() {
        println!("{p}\t{c}\t{d:.4}");
    }


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

