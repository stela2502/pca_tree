//! PCA-Tree: DNA sequence PCA + MST builder

mod encoder;
mod pca;
mod tree;

pub use encoder::OneHotEncoder;
pub use pca::PcaModel;
pub use tree::MstTree;

use ndarray::Array2;
use std::error::Error;

/// Combined PCA + MST pipeline structure.
pub struct PcaTree {
    //encoder: OneHotEncoder,
    pub pca: PcaModel,
    pub tree: MstTree,
}



impl PcaTree {

    /// Build PCA + MST from raw sequences.
    pub fn new(seqs: Vec<String>, k: usize) -> Result<Self, Box<dyn Error>> {

        // Encode sequences numerically
        let encoder = OneHotEncoder::new();
        let encoded = encoder.encode_relative(&seqs)?;

        // Fit PCA
        let mut pca = PcaModel::new(k);
        pca.fit_transform(&encoded)?;

        // Build tree in PCA space
        let tree = MstTree::build(pca.coords());

        Ok(Self {
            //encoder,
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