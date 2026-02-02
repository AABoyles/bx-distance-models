/// Basic distance metrics
use wasm_bindgen::prelude::*;
use crate::utils::*;

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

/// Compute p-distance (simple proportion of differences).
/// Similar to Hamming but with explicit gap handling.
/// Returns the proportion of differing positions where both bases are valid.
#[wasm_bindgen]
pub fn p_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 1.0;
    }

    let mut differences = 0u32;
    let mut compared = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        // Only count positions where both are valid bases
        if !is_valid_base(a) || !is_valid_base(b) {
            continue;
        }

        compared += 1;
        if a != b {
            differences += 1;
        }
    }

    if compared == 0 {
        return 1.0;
    }

    differences as f64 / compared as f64
}
