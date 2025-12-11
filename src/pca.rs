use ndarray::{Array1, Array2, Axis};
use ndarray_linalg::eigh::Eigh;
use std::error::Error;
use ndarray_linalg::UPLO;
#[cfg(feature = "plot")]
use plotters::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter};
use std::fs::{File};
use std::io::Write;
use std::path::Path;


pub struct PcaModel {
    pub k: usize,
    pub mean: Array1<f32>,
    pub components: Array2<f32>,
    pub coords: Array2<f32>,
}

impl PcaModel {

    pub fn new(k: usize) -> Self {
        Self {
            k,
            mean: Array1::zeros(0),
            components: Array2::zeros((0, 0)),
            coords: Array2::zeros((0, 0)),
        }
    }

    /// Write PCA coordinates to TSV (n rows × k columns).
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        self.to_delimited( path, '\t' )
    }

    /// Optional: allow custom separators
    pub fn to_delimited<P: AsRef<Path>>(&self, path: P, sep: char) -> std::io::Result<()> {
        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        for row in self.coords.outer_iter() {
            let mut first = true;
            for v in row {
                if !first {
                    write!(w, "{}", sep)?;
                }
                write!(w, "{:.6}", v)?;
                first = false;
            }
            writeln!(w)?;
        }
        Ok(())
    }
    

    pub fn fit_transform(&mut self, x: &Array2<f32>) -> Result<(), Box<dyn Error>> {
        let (n, p) = x.dim();

        let mean = x.mean_axis(Axis(0)).unwrap();
        let mut centered = x.clone();

        for mut row in centered.outer_iter_mut() {
            row -= &mean;
        }

        let cov = centered.t().dot(&centered) / (n as f32 - 1.0);
        let (eigvals, eigvecs) = cov.eigh(UPLO::Upper)?;

        let mut idx: Vec<_> = (0..eigvals.len()).collect();
        idx.sort_by(|a, b| eigvals[*b].partial_cmp(&eigvals[*a]).unwrap());

        let comps = Array2::from_shape_fn((p, self.k), |(i, j)| eigvecs[(i, idx[j])]);
        let proj = centered.dot(&comps);

        self.mean = mean;
        self.components = comps;
        self.coords = proj;

        Ok(())
    }
    pub fn coords(&self) -> &Array2<f32> {
        &self.coords
    }

    pub fn components(&self) -> &Array2<f32> {
        &self.components
    }
    #[cfg(feature = "plot")]
    pub fn plot_2d_clusters(
        &self,
        tree: &crate::MstTree,
        outfile: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {

        println!("Plotting PCA with cluster coloring → {}", outfile);

        let coords = &self.coords;
        let n = coords.nrows();

        // --- run elbow clustering ---
        let clusters = tree.clusters_elbow(n);

        // map node -> cluster id
        let mut belong = HashMap::<usize, usize>::new();
        for (cid, cluster) in clusters.iter().enumerate() {
            for &node in cluster {
                belong.insert(node, cid);
            }
        }

        // --- build color palette ---
        let palette = [
            RED, BLUE, GREEN, MAGENTA, CYAN, YELLOW,
            RGBColor(255,165,0),   // orange
            RGBColor(128,0,128),   // purple
            RGBColor(0,128,128),   // teal
        ];

        // --- chart setup ---
        let root = BitMapBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = coords.column(0);
        let y = coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA (cluster-colored)", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        // --- draw points ---
        for i in 0..n {
            let color = if let Some(&cid) = belong.get(&i) {
                palette[cid % palette.len()]
            } else {
                RGBColor(160, 160, 160)   // external / orphan
            };

            chart.draw_series([Circle::new(
                (coords[(i,0)], coords[(i,1)]),
                4,
                color.filled()
            )])?;
        }

        root.present()?;
        Ok(())
    }
    #[cfg(feature = "plot")]
    pub fn plot_2d(&self, outfile: &str) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;

        let root = BitMapBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = self.coords.column(0);
        let y = self.coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA projection", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        chart.draw_series(
            x.iter().zip(y.iter())
                .map(|(&x, &y)| Circle::new((x, y), 3, BLUE.filled()))
        )?;

        Ok(())
    }
}

