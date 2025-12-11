use ndarray::{Array2, ArrayView1};
#[allow(dead_code, unused)] // creates a warning otherwise
#[cfg(feature = "plot")]
use plotters::prelude::*;
use std::collections::{VecDeque};

pub struct MstTree {
    pub edges: Vec<(usize, usize, f32)>,
}

impl MstTree {

    pub fn len(&self) -> usize {
        self.edges.len()
    }

    pub fn total_length(&self) -> f32 {
        self.edges.iter().map(|(_, _, d)| d).sum()
    }

    pub fn build(coords: &Array2<f32>) -> Self {
        let n = coords.nrows();
        let mut in_tree = vec![false; n];
        let mut dist = vec![f32::INFINITY; n];
        let mut parent = vec![None; n];

        in_tree[0] = true;

        for i in 1..n {
            dist[i] = Self::euclidean(coords.row(0), coords.row(i));
            parent[i] = Some(0);
        }

        for _ in 1..n - 1 {
            let mut best = None;
            let mut best_d = f32::INFINITY;

            for i in 0..n {
                if !in_tree[i] && dist[i] < best_d {
                    best = Some(i);
                    best_d = dist[i];
                }
            }

            let v = best.unwrap();
            in_tree[v] = true;

            for u in 0..n {
                if !in_tree[u] {
                    let d = Self::euclidean(coords.row(v), coords.row(u));
                    if d < dist[u] {
                        dist[u] = d;
                        parent[u] = Some(v);
                    }
                }
            }
        }

        let mut edges = Vec::new();
        for i in 1..n {
            edges.push((parent[i].unwrap(), i, dist[i]));
        }

        Self {
            edges,
        }
    }


    pub fn clusters_elbow(&self, n_nodes: usize) -> Vec<Vec<usize>> {

        let Some(threshold) = self.elbow_threshold() else {
            return vec![];
        };

        self.clusters_with_cut(n_nodes, threshold)
    }

    pub fn clusters_robust(&self, n_nodes: usize) -> Vec<Vec<usize>> {

        let Some(threshold) = self.robust_threshold_auto() else {
            return vec![];
        };

        self.clusters_with_cut(n_nodes, threshold)
    }

    pub fn clusters_with_cut(&self, n_nodes: usize, max_len: f32) -> Vec<Vec<usize>> {

        let mut adj = vec![Vec::new(); n_nodes];

        for (a, b, d) in &self.edges {
            if *d <= max_len {
                adj[*a].push(*b);
                adj[*b].push(*a);
            }
        }

        let mut visited = vec![false; n_nodes];
        let mut out = Vec::new();

        for i in 0..n_nodes {
            if visited[i] { continue; }

            let mut stack = VecDeque::new();
            let mut comp = Vec::new();

            stack.push_back(i);
            visited[i] = true;

            while let Some(u) = stack.pop_front() {
                comp.push(u);
                for &v in &adj[u] {
                    if !visited[v] {
                        visited[v] = true;
                        stack.push_back(v);
                    }
                }
            }

            out.push(comp);
        }

        out
    }
    /// Automatically chooses clustering threshold using elbow detection.
    pub fn elbow_threshold(&self) -> Option<f32> {

        if self.edges.len() < 2 {
            return None;
        }

        let mut lens: Vec<f32> = self.edges.iter().map(|(_,_,d)| *d).collect();
        lens.sort_by(|a,b| a.partial_cmp(b).unwrap());

        let mut best_i = 0;
        let mut best_gap = 0.0;

        for i in 0..lens.len() - 1 {
            let gap = lens[i+1] - lens[i];
            if gap > best_gap {
                best_gap = gap;
                best_i = i;
            }
        }

        Some(lens[best_i])
    }

    pub fn robust_threshold(&self, k: f32) -> Option<f32> {
        if self.edges.len() < 2 {
            return None;
        }

        let mut x: Vec<f32> = self.edges.iter().map(|(_, _, d)| *d).collect();
        x.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let med = median(&x);

        let dev: Vec<f32> = x.iter().map(|v| (v - med).abs()).collect();
        let mad = median(&dev);

        Some(med + k * mad)
    }

    pub fn robust_threshold_auto(&self) -> Option<f32> {
        if self.edges.len() < 4 {
            return None;
        }

        let mut x: Vec<f32> = self.edges.iter().map(|(_, _, d)| *d).collect();
        x.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // --- median ---
        let med = median(&x);

        // --- MAD ---
        let dev: Vec<f32> = x.iter().map(|v| (v - med).abs()).collect();
        let mad = median(&dev).max(1e-9);

        // --- normalized tail weights ---
        // z-score-like: (x - median) / MAD
        let z: Vec<f32> = x.iter().map(|v| (v - med) / mad).collect();

        // --- detect first big tail rise ---
        // find first value beyond a natural outlier region
        let mut cut = None;

        for i in 0..z.len() {
            // "unlikely under normal" threshold
            if z[i] > 3.5 && i > x.len() / 2 {
                cut = Some(x[i]);
                break;
            }
        }

        // --- fallback: percentile based ---
        if cut.is_none() {
            let idx = ((x.len() as f32) * 0.85) as usize;
            cut = Some(x[idx.min(x.len() - 1)]);
        }

        cut
    }



    pub fn cut(&self, max_len: f32) -> Vec<(usize, usize)> {
        self.edges
            .iter()
            .filter(|(_, _, d)| *d <= max_len)
            .map(|(a, b, _)| (*a, *b))
            .collect()
    }
    #[cfg(feature = "plot")]
    pub fn plot_2d(&self, coords: &ndarray::Array2<f32>, outfile: &str)
        -> Result<(), Box<dyn std::error::Error>>
    {
        use plotters::prelude::*;

        let root = BitMapBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = coords.column(0);
        let y = coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA Tree", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        // Draw edges (lines)
        for &(a, b, _) in &self.edges {
            let pa = (coords[(a, 0)], coords[(a, 1)]);
            let pb = (coords[(b, 0)], coords[(b, 1)]);
            chart.draw_series([PathElement::new(vec![pa, pb], &BLACK)])?;
        }

        // Draw nodes
        chart.draw_series(
            x.iter().zip(y.iter())
                .map(|(&x, &y)| Circle::new((x, y), 3, RED.filled()))
        )?;

        Ok(())
    }

    fn euclidean(a: ArrayView1<f32>, b: ArrayView1<f32>) -> f32 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f32>()
            .sqrt()
    }
}

fn median(v: &[f32]) -> f32 {
    let m = v.len() / 2;
    if v.len() % 2 == 0 {
        (v[m - 1] + v[m]) / 2.0
    } else {
        v[m]
    }
}