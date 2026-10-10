/// Basic distance metrics
use wasm_bindgen::prelude::*;
use crate::counts::Counts;

/// Count the number of nucleotide mismatches (SNPs) between two aligned sequences.
/// Only sites where both sequences hold A, C, G or T (any case, U = T) are compared.
/// Returns the absolute count of differences.
#[wasm_bindgen]
pub fn snps_distance(seq1: &str, seq2: &str) -> f64 {
    snps(&Counts::from_strs(seq1, seq2))
}

/// Compute Hamming distance between two aligned sequences.
/// Returns the proportion of differing positions (0 to 1) among compared sites,
/// or 1.0 if there are none.
#[wasm_bindgen]
pub fn hamming_distance(seq1: &str, seq2: &str) -> f64 {
    hamming(&Counts::from_strs(seq1, seq2))
}

/// Compute p-distance (simple proportion of differences).
/// Identical to `hamming_distance`: both compare only sites where both sequences hold a base.
#[wasm_bindgen]
pub fn p_distance(seq1: &str, seq2: &str) -> f64 {
    hamming(&Counts::from_strs(seq1, seq2))
}

pub(crate) fn snps(c: &Counts) -> f64 {
    c.mismatches()
}

pub(crate) fn hamming(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        1.0
    } else {
        c.mismatches() / c.sites
    }
}
