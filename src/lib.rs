//! ClonoMap: scalable inference of clonal structure in large AIRR-seq clones.
//!
//! ClonoMap infers subclonal structure within B-cell receptor (BCR) clones
//! using mutation-aware geometric clustering rather than full phylogenetic
//! reconstruction. It is designed to scale to clones containing tens of
//! thousands of sequences and integrates seamlessly with Change-O outputs.

mod encoder;
mod pca;
mod tree;
mod clone_data;

pub use encoder::OneHotEncoder;
pub use pca::PcaModel;
pub use tree::MstTree;
pub use clone_data::CloneData;

use ndarray::Array2;
use std::error::Error;
use std::collections::HashMap;


/// Combined PCA + MST pipeline structure.
pub struct PcaTree {
    pub encoder: OneHotEncoder,
    pub pca: PcaModel,
    pub tree: MstTree,
}



impl PcaTree {

    /// Build PCA + MST from raw sequences.
    pub fn new(seqs: Vec<String>, k: usize) -> Result<Self, Box<dyn Error>> {

        /*let mut map: HashMap<String, usize> = HashMap::new();

        for s in seqs {
            *map.entry(s).or_insert(0) += 1;
        }

        // Unique sequences
        let unique: Vec<String> = map.keys().cloned().collect();
        // Counts for each unique seq
        let counts: Vec<usize> = unique.iter().map(|s| map[s]).collect();
        */

        // Encode sequences numerically
        let mut encoder = OneHotEncoder::new();
        println!("PcaTree::new - I got {} sequences", seqs.len() );

        let encoded = encoder.encode_relative(&seqs)?;

        println!("      and I encoded {} unique amino acid respresentations ({} cols)", encoded.nrows(),  encoded.ncols() );

        // Fit PCA
        let mut pca = PcaModel::new(k);
        pca.fit_transform(&encoded)?;

        // Build tree in PCA space
        let tree = MstTree::build(&pca);

        Ok(Self {
            encoder,
            pca,
            tree,
        })
    }

    /// PCA coordinates accessor
    pub fn coords(&self) -> &Array2<f32> {
        self.pca.coords()
    }

    /// Tree edge list
    pub fn tree(&self) -> &Vec<(usize, usize, f32)> {
        &self.tree.edges

    }
}
