use ndarray::{Array2, ArrayView1};
#[allow(dead_code, unused)] // creates a warning otherwise
#[cfg(feature = "plot")]
use plotters::prelude::*;
use std::collections::{VecDeque};
use std::io::{BufWriter};
use std::fs::{File};
use std::io::Write;
use std::path::Path;
use rayon::prelude::*;
use crate::CloneData;

use crate::PcaModel;

pub struct MstTree {
    pub edges: Vec<(usize, usize, f32)>,
}

impl MstTree {

    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Write MST edges as TSV: parent<TAB>child<TAB>distance
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        self.to_delimited(path, '\t')
    }

    /// Write MST edges using a custom delimiter.
    pub fn to_delimited<P: AsRef<Path>>(&self, path: P, sep: char) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};

        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        for (p, c, d) in &self.edges {
            writeln!(w, "{}{}{}{}.{:.6}", p, sep, c, sep, d)?;
        }

        Ok(())
    }

    pub fn total_length(&self) -> f32 {
        self.edges.iter().map(|(_, _, d)| d).sum()
    }

    pub fn build(model: &PcaModel) -> Self {
        use rayon::prelude::*;
        let coords = &model.coords;

        let n = coords.nrows();
        if n == 0 {
            return Self { edges: Vec::new() };
        }

        // MST state
        let mut in_tree = vec![false; n];
        let mut dist = vec![f32::INFINITY; n];
        let mut parent = vec![None; n];

        in_tree[0] = true;

        // Initialize distances relative to node 0
        for i in 1..n {
            dist[i] = Self::euclidean(coords.row(0), coords.row(i));
            parent[i] = Some(0);
        }

        // --- main MST loop ---
        for _ in 1..(n - 1) {
            // Find best vertex outside MST (sequential O(n))
            let mut best = None;
            let mut best_d = f32::INFINITY;

            for i in 0..n {
                if !in_tree[i] && dist[i] < best_d {
                    best = Some(i);
                    best_d = dist[i];
                }
            }

            let v = best.expect("MST cannot proceed – no available node");
            in_tree[v] = true;

            let v_row = coords.row(v);

            // --- PARALLEL relax edges ---
            // Instead of updating dist[u] directly (data race)
            // we compute all improvements first, then apply them.
            let updates: Vec<(usize, f32, usize)> = 
                (0..n)
                    .into_par_iter()
                    .filter_map(|u| {
                        if in_tree[u] {
                            return None;
                        }

                        let d = Self::euclidean(v_row, coords.row(u));
                        if d < dist[u] {
                            Some((u, d, v)) // (node, new-dist, parent)
                        } else {
                            None
                        }
                    })
                    .collect();

            // Apply computed updates (sequential, very cheap)
            for (u, newd, pv) in updates {
                dist[u] = newd;
                parent[u] = Some(pv);
            }
        }

        // Collect MST edges
        let mut edges = Vec::with_capacity(n - 1);
        for i in 1..n {
            edges.push((parent[i].unwrap(), i, dist[i]));
        }

        Self { edges }
    }

    /// Find sparse root: the node located in the lowest-density PCA region.
    /// Uses the distance to the k-th nearest neighbor.
    pub fn find_sparse_root(model: &PcaModel, k: usize) -> usize {
        let coords = &model.coords;
        let n = coords.nrows();
        let mut scores = vec![0.0; n];

        for i in 0..n {
            // Distance to all other nodes
            let mut dists: Vec<f32> = (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    coords.row(i)
                        .iter()
                        .zip(coords.row(j).iter())
                        .map(|(x, y)| (x - y).powi(2))
                        .sum::<f32>()
                        .sqrt()
                })
                .collect();

            // Sort and store k-th smallest distance
            dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
            scores[i] = dists[k.min(dists.len() - 1)];
        }

        // Node with the largest kNN radius = sparsest region = predicted root
        scores
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0
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

    /// Re-root the MST at the given node index
    pub fn reroot(&self, n: usize, root: usize) -> MstTree {
        let mut adj = vec![Vec::new(); n];

        // Build adjacency
        for &(p, c, d) in &self.edges {
            adj[p].push((c, d));
            adj[c].push((p, d));
        }

        let mut parent = vec![None; n];
        let mut dist_to_parent = vec![0.0f32; n];

        let mut queue = VecDeque::new();
        queue.push_back(root);
        parent[root] = Some(root); // mark root

        // BFS to orient edges
        while let Some(v) = queue.pop_front() {
            for &(nbr, d) in &adj[v] {
                if parent[nbr].is_none() {
                    parent[nbr] = Some(v);
                    dist_to_parent[nbr] = d;
                    queue.push_back(nbr);
                }
            }
        }

        // Rebuild edges (skip the root)
        let mut new_edges = Vec::new();
        for i in 0..n {
            if i == root {
                continue;
            }
            let p = parent[i].unwrap();
            new_edges.push((p, i, dist_to_parent[i]));
        }

        MstTree { edges: new_edges }
    }

    pub fn to_newick(&self, n: usize, root: usize, labels: &CloneData ) -> String {
        let mut children = vec![Vec::new(); n];
        for &(p, c, d) in &self.edges {
            children[p].push((c, d));
        }

        fn build(idx: usize, children: &Vec<Vec<(usize, f32)>>, aa_labels:&[String] ) -> String {
            if children[idx].is_empty() {
                format!("{}", aa_labels[idx])
            } else {
                let inner: Vec<String> = children[idx]
                    .iter()
                    .map(|(c, d)| format!("{}:{:.4}", build(*c, children, aa_labels), d))
                    .collect();

                format!("({})", inner.join(","))
            }
        }
        let keys = labels.aa_with_count_labels();

        format!("{};", build(root, &children, &keys))
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