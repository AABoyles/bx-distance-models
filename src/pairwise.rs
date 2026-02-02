/// Pairwise distance computation functions
use wasm_bindgen::prelude::*;
use js_sys::{Array, Float64Array, Object, Reflect};
use crate::metrics::{hamming_distance, tn93_distance, snps_distance};

/// Compute all pairwise distances for a list of sequences.
/// Takes an array of {id: string, seq: string} objects and returns
/// an array of {source: string, target: string, distance: number} objects.
///
/// This is the main entry point for bulk distance computation.
#[wasm_bindgen]
pub fn compute_pairwise_distances(nodes: &Array, metric: &str) -> Array {
    let n = nodes.length() as usize;
    let result = Array::new();

    // Extract sequences into Vec for faster access
    let mut seqs: Vec<(String, String)> = Vec::with_capacity(n);

    for i in 0..n {
        let node = nodes.get(i as u32);
        if let (Ok(id), Ok(seq)) = (
            Reflect::get(&node, &JsValue::from_str("id")),
            Reflect::get(&node, &JsValue::from_str("seq")),
        ) {
            if let (Some(id_str), Some(seq_str)) = (id.as_string(), seq.as_string()) {
                seqs.push((id_str, seq_str));
            }
        }
    }

    let use_tn93 = metric == "tn93";
    let seq_count = seqs.len();

    for i in 0..seq_count {
        for j in (i + 1)..seq_count {
            let dist = if use_tn93 {
                tn93_distance(&seqs[i].1, &seqs[j].1)
            } else {
                hamming_distance(&seqs[i].1, &seqs[j].1)
            };

            let edge = Object::new();
            let _ = Reflect::set(&edge, &JsValue::from_str("source"), &JsValue::from_str(&seqs[i].0));
            let _ = Reflect::set(&edge, &JsValue::from_str("target"), &JsValue::from_str(&seqs[j].0));
            let _ = Reflect::set(&edge, &JsValue::from_str("distance"), &JsValue::from_f64(dist));
            result.push(&edge);
        }
    }

    result
}

/// Compute pairwise distances using flat arrays for maximum performance.
/// Takes:
/// - _ids: Array of node ID strings (unused, provided for consistency)
/// - seqs: Array of sequence strings
/// - metric: "hamming" or "tn93"
///
/// Returns a Float64Array with distances in row-major upper triangular order.
/// For n nodes, returns n*(n-1)/2 distances.
#[wasm_bindgen]
pub fn compute_pairwise_distances_flat(_ids: &Array, seqs: &Array, metric: &str) -> Float64Array {
    let n = seqs.length() as usize;
    let num_pairs = n * (n - 1) / 2;

    // Extract sequences into Vec<String> for faster access
    let seq_vec: Vec<String> = (0..n)
        .filter_map(|i| seqs.get(i as u32).as_string())
        .collect();

    if seq_vec.len() != n {
        return Float64Array::new_with_length(0);
    }

    let use_tn93 = metric == "tn93";
    let mut distances = vec![0.0f64; num_pairs];
    let mut idx = 0;

    for i in 0..n {
        for j in (i + 1)..n {
            distances[idx] = if use_tn93 {
                tn93_distance(&seq_vec[i], &seq_vec[j])
            } else {
                hamming_distance(&seq_vec[i], &seq_vec[j])
            };
            idx += 1;
        }
    }

    Float64Array::from(distances.as_slice())
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
/// This efficiently computes all metrics in a single pass.
#[wasm_bindgen]
pub fn compute_pairwise_distances_flat_all(_ids: &Array, seqs: &Array) -> Float64Array {
    let n = seqs.length() as usize;
    let num_pairs = n * (n - 1) / 2;

    // Extract sequences into Vec<String> for faster access
    let seq_vec: Vec<String> = (0..n)
        .filter_map(|i| seqs.get(i as u32).as_string())
        .collect();

    if seq_vec.len() != n {
        return Float64Array::new_with_length(0);
    }

    // Allocate space for all three metrics
    let total_len = num_pairs * 3;
    let mut distances = vec![0.0f64; total_len];
    let mut idx = 0;

    for i in 0..n {
        for j in (i + 1)..n {
            // Compute all metrics in one pass
            let seq1 = &seq_vec[i];
            let seq2 = &seq_vec[j];

            distances[idx] = hamming_distance(seq1, seq2);
            distances[idx + num_pairs] = tn93_distance(seq1, seq2);
            distances[idx + num_pairs * 2] = snps_distance(seq1, seq2);

            idx += 1;
        }
    }

    Float64Array::from(distances.as_slice())
}

/// Get the index into the flat distance array for a given (i, j) pair.
/// Assumes i < j and n is the total number of nodes.
#[wasm_bindgen]
pub fn get_distance_index(i: u32, j: u32, n: u32) -> u32 {
    // Upper triangular index formula
    let i = i as usize;
    let j = j as usize;
    let n = n as usize;
    (i * (2 * n - i - 1) / 2 + (j - i - 1)) as u32
}
