//! Batch API: encode a set of aligned sequences once, then compute any metric over all pairs.
use std::ops::Range;

use js_sys::{Float64Array, Uint32Array};
use wasm_bindgen::prelude::*;

use crate::counts::{count_pair, Counts, Resolutions};
use crate::encoding::{Ambiguity, Encoded};
use crate::metrics::{Metric, METRIC_NAMES};

/// A set of aligned nucleotide sequences, encoded once for fast pairwise comparison.
///
/// Pairs are ordered like the upper triangle of the distance matrix, row by row:
/// (0,1), (0,2), …, (0,n-1), (1,2), …. Every method can be limited to a range of rows, so the
/// work can be split into chunks (progress, cancellation) or across Web Workers.
#[wasm_bindgen]
pub struct SequenceSet {
    seqs: Vec<Encoded>,
    /// Ambiguity resolution table; `None` in "skip" mode.
    resolutions: Option<Resolutions>,
}

/// Pairs within a distance threshold, as parallel arrays: `source[k] < target[k]` are sequence
/// indices and `distance[k]` their distance. Each getter returns a fresh copy; read it once.
#[wasm_bindgen]
#[derive(Debug, Default, PartialEq)]
pub struct Links {
    pub(crate) source: Vec<u32>,
    pub(crate) target: Vec<u32>,
    pub(crate) distance: Vec<f64>,
}

#[wasm_bindgen]
impl Links {
    #[wasm_bindgen(getter)]
    pub fn length(&self) -> u32 {
        self.distance.len() as u32
    }

    #[wasm_bindgen(getter)]
    pub fn source(&self) -> Uint32Array {
        Uint32Array::from(&self.source[..])
    }

    #[wasm_bindgen(getter)]
    pub fn target(&self) -> Uint32Array {
        Uint32Array::from(&self.target[..])
    }

    #[wasm_bindgen(getter)]
    pub fn distance(&self) -> Float64Array {
        Float64Array::from(&self.distance[..])
    }
}

#[wasm_bindgen]
impl SequenceSet {
    /// Encodes `seqs` once. `ambiguity` decides how sites holding IUPAC ambiguity codes
    /// (R, Y, S, W, K, M, B, D, H, V) are compared:
    /// - `"skip"` (default): ignore them. Every pair uses the fast bit-plane kernel.
    /// - `"resolve"`: count a site as a match if the two codes share a base, otherwise average.
    /// - `"average"`: average over every pair of bases the two codes allow.
    ///
    /// With `"resolve"` or `"average"`, only pairs involving a sequence that actually contains
    /// ambiguity codes pay for the slower fallback. `N`, gaps and other characters are always
    /// skipped.
    #[wasm_bindgen(constructor)]
    pub fn new(seqs: Vec<String>, ambiguity: Option<String>) -> Result<SequenceSet, JsError> {
        let mode = match ambiguity.as_deref() {
            None => Ambiguity::Skip,
            Some(name) => Ambiguity::parse(name).ok_or_else(|| {
                JsError::new(&format!("unknown ambiguity mode '{name}' (expected skip, resolve or average)"))
            })?,
        };
        Ok(SequenceSet::from_strs(&seqs, mode))
    }

    /// Number of sequences.
    #[wasm_bindgen(getter)]
    pub fn length(&self) -> u32 {
        self.seqs.len() as u32
    }

    /// Number of sequences that contain ambiguity codes and use the fallback kernel
    /// (always 0 in `"skip"` mode).
    pub fn ambiguous_count(&self) -> u32 {
        self.seqs.iter().filter(|s| s.ambiguity.is_some()).count() as u32
    }

    /// Distance between sequences `i` and `j`.
    pub fn distance(&self, metric: &str, i: u32, j: u32) -> Result<f64, JsError> {
        let metric = parse_metric(metric).map_err(|e| JsError::new(&e))?;
        let n = self.len();
        if i as usize >= n || j as usize >= n {
            return Err(JsError::new(&format!("index out of range for {n} sequences")));
        }
        Ok(self.pair(metric, i as usize, j as usize))
    }

    /// Dense upper-triangle distances for rows `[row_start, row_end)` (default: all rows),
    /// in pair order. Concatenating the results for consecutive row ranges gives the full
    /// triangle; `get_distance_index(i, j, n)` maps a pair to its position in it.
    pub fn distances(&self, metric: &str, row_start: Option<u32>, row_end: Option<u32>) -> Result<Float64Array, JsError> {
        let metric = parse_metric(metric).map_err(|e| JsError::new(&e))?;
        let rows = self.row_range(row_start, row_end).map_err(|e| JsError::new(&e))?;
        self.dense_js(&[metric], rows).map_err(|e| JsError::new(&e))
    }

    /// Pairs in rows `[row_start, row_end)` whose distance is `<= threshold`.
    ///
    /// For `snps`, `hamming`, `p_distance`, `jc69`, `k2p` and `tn93`, a pair stops being scanned
    /// as soon as it can no longer pass, which makes sparse output much faster than
    /// `distances()` when most pairs are far apart.
    pub fn links(&self, metric: &str, threshold: f64, row_start: Option<u32>, row_end: Option<u32>) -> Result<Links, JsError> {
        let metric = parse_metric(metric).map_err(|e| JsError::new(&e))?;
        let rows = self.row_range(row_start, row_end).map_err(|e| JsError::new(&e))?;
        Ok(self.sparse(metric, threshold, rows))
    }

    /// `parts + 1` row boundaries that split the pairs into `parts` chunks of nearly equal
    /// size: chunk `k` covers rows `[splits[k], splits[k + 1])`.
    pub fn row_splits(&self, parts: u32) -> Uint32Array {
        Uint32Array::from(&row_splits(self.len(), parts as usize)[..])
    }
}

impl SequenceSet {
    pub(crate) fn from_strs<S: AsRef<str>>(seqs: &[S], mode: Ambiguity) -> SequenceSet {
        let resolutions = (mode != Ambiguity::Skip).then(|| Resolutions::new(mode));
        let track = resolutions.is_some();
        let seqs = seqs.iter().map(|s| Encoded::new(s.as_ref().as_bytes(), track)).collect();
        SequenceSet { seqs, resolutions }
    }

    pub(crate) fn len(&self) -> usize {
        self.seqs.len()
    }

    #[inline]
    pub(crate) fn counts(&self, i: usize, j: usize, max_mismatches: f64) -> Option<Counts> {
        count_pair(&self.seqs[i], &self.seqs[j], self.resolutions.as_ref(), max_mismatches)
    }

    pub(crate) fn pair(&self, metric: Metric, i: usize, j: usize) -> f64 {
        let c = self.counts(i, j, f64::INFINITY).expect("an unbounded scan always completes");
        metric.eval(&c)
    }

    fn row_range(&self, start: Option<u32>, end: Option<u32>) -> Result<Range<usize>, String> {
        let n = self.len();
        let (start, end) = (start.map_or(0, |s| s as usize), end.map_or(n, |e| e as usize));
        if start > end || end > n {
            return Err(format!("invalid row range {start}..{end} for {n} sequences"));
        }
        Ok(start..end)
    }

    /// Computes `metrics` for every pair in `rows`, row by row. `sink(offset, values)` receives
    /// each row's values for one metric; `offset` indexes a metric-major output array
    /// (all pairs of metric 0, then all pairs of metric 1, ...).
    fn dense_rows(&self, metrics: &[Metric], rows: Range<usize>, mut sink: impl FnMut(usize, &[f64])) {
        let n = self.len();
        let total = pairs_in_rows(n, rows.clone()) as usize;
        let mut row = Vec::new();
        let mut offset = 0;
        for i in rows {
            let width = n - 1 - i;
            row.resize(width * metrics.len(), 0.0);
            for j in i + 1..n {
                let c = self.counts(i, j, f64::INFINITY).expect("an unbounded scan always completes");
                for (m, metric) in metrics.iter().enumerate() {
                    row[m * width + (j - i - 1)] = metric.eval(&c);
                }
            }
            for m in 0..metrics.len() {
                sink(m * total + offset, &row[m * width..(m + 1) * width]);
            }
            offset += width;
        }
    }

    /// Dense output as a `Vec` (metric-major).
    #[cfg(test)]
    pub(crate) fn dense(&self, metrics: &[Metric], rows: Range<usize>) -> Vec<f64> {
        let total = pairs_in_rows(self.len(), rows.clone()) as usize;
        let mut out = vec![0.0; total * metrics.len()];
        self.dense_rows(metrics, rows, |at, values| out[at..at + values.len()].copy_from_slice(values));
        out
    }

    /// Dense output written straight into a JS `Float64Array`, one row at a time, so peak
    /// memory is the result plus one row (not the result twice).
    pub(crate) fn dense_js(&self, metrics: &[Metric], rows: Range<usize>) -> Result<Float64Array, String> {
        let len = pairs_in_rows(self.len(), rows.clone()) * metrics.len() as u64;
        if len * 8 > u32::MAX as u64 {
            return Err(format!(
                "{len} distances do not fit in WebAssembly memory; use links() with a threshold or smaller row ranges"
            ));
        }
        let out = Float64Array::new_with_length(len as u32);
        self.dense_rows(metrics, rows, |at, values| {
            out.subarray(at as u32, (at + values.len()) as u32).copy_from(values)
        });
        Ok(out)
    }

    pub(crate) fn sparse(&self, metric: Metric, threshold: f64, rows: Range<usize>) -> Links {
        let mut links = Links::default();
        for i in rows {
            let max_i = self.seqs[i].informative;
            for j in i + 1..self.len() {
                let bound = metric.mismatch_bound(threshold, max_i.min(self.seqs[j].informative));
                if let Some(c) = self.counts(i, j, bound) {
                    let d = metric.eval(&c);
                    if d <= threshold {
                        links.source.push(i as u32);
                        links.target.push(j as u32);
                        links.distance.push(d);
                    }
                }
            }
        }
        links
    }
}

fn parse_metric(name: &str) -> Result<Metric, String> {
    Metric::parse(name).ok_or_else(|| format!("unknown metric '{name}' (expected one of: {METRIC_NAMES})"))
}

/// Number of pairs before row `i` of the upper triangle of an `n`-sequence matrix.
pub(crate) fn pairs_before_row(n: usize, i: usize) -> u64 {
    let (n, i) = (n as u64, i as u64);
    if i == 0 { 0 } else { i * (2 * n - i - 1) / 2 }
}

pub(crate) fn pairs_in_rows(n: usize, rows: Range<usize>) -> u64 {
    pairs_before_row(n, rows.end) - pairs_before_row(n, rows.start)
}

/// `parts + 1` row boundaries splitting the triangle into chunks with nearly equal pair counts.
pub(crate) fn row_splits(n: usize, parts: usize) -> Vec<u32> {
    let parts = parts.max(1);
    let total = pairs_before_row(n, n);
    let mut splits = Vec::with_capacity(parts + 1);
    splits.push(0);
    let (mut row, mut done) = (0usize, 0u64);
    for k in 1..parts {
        let target = total * k as u64 / parts as u64;
        while row < n && done < target {
            done += (n - 1 - row) as u64;
            row += 1;
        }
        splits.push(row as u32);
    }
    splits.push(n as u32);
    splits
}
