use ndarray::Array2;
use std::error::Error;

pub struct OneHotEncoder;

impl OneHotEncoder {
    pub fn new() -> Self {
        Self
    }

    pub fn encode_batch(&self, sequences: &[String]) -> Result<Array2<f32>, Box<dyn Error>> {
        if sequences.is_empty() {
            return Err("No sequences provided".into());
        }

        let len = sequences[0].len();

        for (i, s) in sequences.iter().enumerate() {
            if s.len() != len {
                return Err(format!(
                    "Sequence length mismatch at index {}: expected {}, got {}",
                    i, len, s.len()
                ).into());
            }
        }

        let n = sequences.len();
        let d = 4 * len;

        let mut x = Array2::<f32>::zeros((n, d));

        for (i, seq) in sequences.iter().enumerate() {
            for (pos, base) in seq.chars().enumerate() {
                let idx = match base {
                    'A' | 'a' => 0,
                    'C' | 'c' => 1,
                    'G' | 'g' => 2,
                    'T' | 't' => 3,
                    _ => return Err(format!("Invalid base {base}").into()),
                };
                x[[i, 4 * pos + idx]] = 1.0;
            }
        }

        Ok(x)
    }

    pub fn encode_relative(&self, sequences: &[String]) -> Result<Array2<f32>, Box<dyn Error>> {
        if sequences.is_empty() {
            return Err("No sequences provided".into());
        }

        let len = sequences[0].len();

        for (i, s) in sequences.iter().enumerate() {
            if s.len() != len {
                return Err(format!(
                    "Sequence length mismatch at index {}: expected {}, got {}",
                    i, len, s.len()
                ).into());
            }
        }

        let n = sequences.len();
        let mut x = Array2::<f32>::zeros((n, len));

        // consensus per column (ignore gaps)
        let mut consensus = Vec::with_capacity(len);

        for col in 0..len {
            let mut counts = [0u32; 4];

            for s in sequences {
                match s.as_bytes()[col] {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' => counts[3] += 1,
                    _ => {}
                }
            }

            let (idx, _) = counts.iter().enumerate().max_by_key(|(_, c)| *c).unwrap();
            consensus.push(b"ACGT"[idx]);
        }

        // Encode
        for (i, seq) in sequences.iter().enumerate() {
            for (j, b) in seq.as_bytes().iter().enumerate() {
                x[[i, j]] = match *b {
                    b'.' | b'-' => -1.0,
                    _ => {
                        if *b == consensus[j] { 0.0 } else { 1.0 }
                    }
                };
            }
        }

        Ok(x)
    }

}

