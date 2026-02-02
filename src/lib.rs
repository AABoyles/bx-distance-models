use wasm_bindgen::prelude::*;
use js_sys::{Array, Float64Array, Object, Reflect};

#[cfg(feature = "console_error_panic_hook")]
use console_error_panic_hook;

/// Initialize panic hook for better error messages in browser console
#[wasm_bindgen(start)]
pub fn init() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

/// Count the number of nucleotide mismatches (SNPs) between two aligned sequences.
/// Ignores positions where either sequence has a gap '-' or 'N'.
/// Returns the absolute count of differences.
#[wasm_bindgen]
pub fn snps_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 0.0;
    }

    let mut mismatches = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        // Skip gaps and N's
        if a == b'-' || b == b'-' || a == b'N' || b == b'N' {
            continue;
        }

        if a != b {
            mismatches += 1;
        }
    }

    mismatches as f64
}

/// Compute Hamming distance between two aligned sequences.
/// Returns the proportion of differing positions (0 to 1).
/// Ignores positions where either sequence has a gap '-' or 'N'.
#[wasm_bindgen]
pub fn hamming_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 1.0;
    }

    let mut mismatches = 0u32;
    let mut compared = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        // Skip gaps and N's
        if a == b'-' || b == b'-' || a == b'N' || b == b'N' {
            continue;
        }

        compared += 1;
        if a != b {
            mismatches += 1;
        }
    }

    if compared == 0 {
        1.0
    } else {
        mismatches as f64 / compared as f64
    }
}

/// Compute Tamura-Nei 93 (TN93) distance between two aligned sequences.
/// This is a more sophisticated nucleotide substitution model that accounts
/// for different rates of transitions and transversions.
#[wasm_bindgen]
pub fn tn93_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 1.0;
    }

    let mut total_sites = 0u32;
    let mut transitions1 = 0u32; // A<->G (purine transitions)
    let mut transitions2 = 0u32; // C<->T (pyrimidine transitions)
    let mut transversions = 0u32;

    let mut freq_a = 0u32;
    let mut freq_c = 0u32;
    let mut freq_g = 0u32;
    let mut freq_t = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        // Only count valid bases
        if !is_valid_base(a) || !is_valid_base(b) {
            continue;
        }

        total_sites += 1;

        // Count frequencies
        match a {
            b'A' => freq_a += 1,
            b'C' => freq_c += 1,
            b'G' => freq_g += 1,
            b'T' => freq_t += 1,
            _ => {}
        }
        match b {
            b'A' => freq_a += 1,
            b'C' => freq_c += 1,
            b'G' => freq_g += 1,
            b'T' => freq_t += 1,
            _ => {}
        }

        if a == b {
            continue;
        }

        // Classify substitution type
        if is_purine(a) && is_purine(b) {
            transitions1 += 1;
        } else if is_pyrimidine(a) && is_pyrimidine(b) {
            transitions2 += 1;
        } else {
            transversions += 1;
        }
    }

    if total_sites == 0 {
        return 1.0;
    }

    // Normalize frequencies
    let total2 = (total_sites * 2) as f64;
    let g_a = freq_a as f64 / total2;
    let g_c = freq_c as f64 / total2;
    let g_g = freq_g as f64 / total2;
    let g_t = freq_t as f64 / total2;

    let g_r = g_a + g_g; // purine frequency
    let g_y = g_c + g_t; // pyrimidine frequency

    if g_r == 0.0 || g_y == 0.0 {
        return hamming_distance(seq1, seq2);
    }

    let p1 = transitions1 as f64 / total_sites as f64;
    let p2 = transitions2 as f64 / total_sites as f64;
    let q = transversions as f64 / total_sites as f64;

    // TN93 formula
    let a1 = (2.0 * g_a * g_g) / g_r;
    let a2 = (2.0 * g_c * g_t) / g_y;
    let b = 2.0 * (g_r * g_y - (g_a * g_g * g_y) / g_r - (g_c * g_t * g_r) / g_y);

    let term1 = 1.0 - p1 / a1 - q / (2.0 * g_r);
    let term2 = 1.0 - p2 / a2 - q / (2.0 * g_y);
    let term3 = 1.0 - q / (2.0 * g_r * g_y);

    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 {
        // Distance is too large for the model; fall back to Hamming
        return hamming_distance(seq1, seq2);
    }

    let d = -a1 * term1.ln() - a2 * term2.ln() - b * term3.ln();

    d.max(0.0).min(1.0)
}

#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

#[inline]
fn is_purine(b: u8) -> bool {
    matches!(b, b'A' | b'G')
}

#[inline]
fn is_pyrimidine(b: u8) -> bool {
    matches!(b, b'C' | b'T')
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_identical() {
        assert_eq!(hamming_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_hamming_different() {
        assert_eq!(hamming_distance("AAAA", "TTTT"), 1.0);
    }

    #[test]
    fn test_hamming_partial() {
        assert_eq!(hamming_distance("AAGT", "ACGT"), 0.25);
    }

    #[test]
    fn test_hamming_with_gaps() {
        // With gap, only 3 positions compared, 1 mismatch
        assert!((hamming_distance("A-GT", "ACGT") - 0.333).abs() < 0.01);
    }

    #[test]
    fn test_tn93_identical() {
        assert_eq!(tn93_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_snps_identical() {
        assert_eq!(snps_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_snps_different() {
        assert_eq!(snps_distance("AAAA", "TTTT"), 4.0);
    }

    #[test]
    fn test_snps_with_gaps() {
        assert_eq!(snps_distance("A-GT", "ACGT"), 1.0);
    }
}
