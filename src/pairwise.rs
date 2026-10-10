/// Pairwise distance computation functions.
///
/// These predate `SequenceSet` and keep their signatures; they now run on the same encoded
/// kernel. Prefer `SequenceSet` for new code: it encodes once, supports every metric, row
/// ranges, ambiguity modes and thresholded output.
use wasm_bindgen::prelude::*;
use js_sys::{Array, Float64Array, Object, Reflect};
use crate::encoding::Ambiguity;
use crate::metrics::Metric;
use crate::sequence_set::{pairs_before_row, SequenceSet};

/// Any metric name `SequenceSet` accepts; unknown names fall back to Hamming as before.
fn legacy_metric(name: &str) -> Metric {
    Metric::parse(name).unwrap_or(Metric::Hamming)
}

/// Extracts every element of `seqs` as a string, or `None` if any is not a string.
fn strings(seqs: &Array) -> Option<Vec<String>> {
    seqs.iter().map(|s| s.as_string()).collect()
}

/// Compute all pairwise distances for a list of sequences.
/// Takes an array of {id: string, seq: string} objects and returns
/// an array of {source: string, target: string, distance: number} objects.
///
/// Creating one JS object per pair costs far more than the distance itself;
/// `SequenceSet.distances()` or `SequenceSet.links()` avoid it.
#[wasm_bindgen]
pub fn compute_pairwise_distances(nodes: &Array, metric: &str) -> Array {
    let (key_id, key_seq) = (JsValue::from_str("id"), JsValue::from_str("seq"));
    let mut ids = Vec::with_capacity(nodes.length() as usize);
    let mut seqs = Vec::with_capacity(nodes.length() as usize);
    for node in nodes.iter() {
        if let (Ok(id), Ok(seq)) = (Reflect::get(&node, &key_id), Reflect::get(&node, &key_seq)) {
            if let (true, Some(seq)) = (id.is_string(), seq.as_string()) {
                ids.push(id);
                seqs.push(seq);
            }
        }
    }

    let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
    let metric = legacy_metric(metric);
    let (key_source, key_target, key_distance) =
        (JsValue::from_str("source"), JsValue::from_str("target"), JsValue::from_str("distance"));
    let result = Array::new();
    for i in 0..ids.len() {
        for j in (i + 1)..ids.len() {
            let edge = Object::new();
            let _ = Reflect::set(&edge, &key_source, &ids[i]);
            let _ = Reflect::set(&edge, &key_target, &ids[j]);
            let _ = Reflect::set(&edge, &key_distance, &JsValue::from_f64(set.pair(metric, i, j)));
            result.push(&edge);
        }
    }
    result
}

/// Compute pairwise distances using flat arrays for maximum performance.
/// Takes:
/// - _ids: Array of node ID strings (unused, provided for consistency)
/// - seqs: Array of sequence strings
/// - metric: any metric name accepted by `SequenceSet` (unknown names use "hamming")
///
/// Returns a Float64Array with distances in row-major upper triangular order.
/// For n nodes, returns n*(n-1)/2 distances.
#[wasm_bindgen]
pub fn compute_pairwise_distances_flat(_ids: &Array, seqs: &Array, metric: &str) -> Float64Array {
    let Some(seqs) = strings(seqs) else {
        return Float64Array::new_with_length(0);
    };
    let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
    set.dense_js(&[legacy_metric(metric)], 0..set.len())
        .unwrap_or_else(|e| wasm_bindgen::throw_str(&e))
}

/// Compute all pairwise distances (hamming, tn93, and snps) using flat arrays.
/// Takes:
/// - _ids: Array of node ID strings (unused, provided for consistency)
/// - seqs: Array of sequence strings
///
/// Returns three Float64Arrays in sequence:
/// - hammingDistances[0..num_pairs)
/// - tn93Distances[num_pairs..2*num_pairs)
/// - snpsDistances[2*num_pairs..3*num_pairs)
///
/// Each pair is scanned once for all three metrics.
#[wasm_bindgen]
pub fn compute_pairwise_distances_flat_all(_ids: &Array, seqs: &Array) -> Float64Array {
    let Some(seqs) = strings(seqs) else {
        return Float64Array::new_with_length(0);
    };
    let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
    set.dense_js(&[Metric::Hamming, Metric::Tn93, Metric::Snps], 0..set.len())
        .unwrap_or_else(|e| wasm_bindgen::throw_str(&e))
}

/// Get the index into the flat distance array for a given (i, j) pair.
/// Assumes i < j and n is the total number of nodes.
#[wasm_bindgen]
pub fn get_distance_index(i: u32, j: u32, n: u32) -> u32 {
    // Upper triangular index formula, in 64 bits so the intermediate product cannot overflow
    (pairs_before_row(n as usize, i as usize) + (j - i - 1) as u64) as u32
}
