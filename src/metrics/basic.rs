/// Basic distance metrics
use wasm_bindgen::prelude::*;
use crate::counts::Mismatches;

/// Count the number of nucleotide mismatches (SNPs) between two aligned sequences.
/// Only sites where both sequences hold A, C, G or T (any case, U = T) are compared.
/// Returns the absolute count of differences.
#[wasm_bindgen]
pub fn snps_distance(seq1: &str, seq2: &str) -> f64 {
    snps(&Mismatches::from_strs(seq1, seq2))
}

/// Compute Hamming distance between two aligned sequences.
/// Returns the proportion of differing positions (0 to 1) among compared sites,
/// or NaN if there are none (the distance is undefined).
#[wasm_bindgen]
pub fn hamming_distance(seq1: &str, seq2: &str) -> f64 {
    hamming(&Mismatches::from_strs(seq1, seq2))
}

/// Compute p-distance (simple proportion of differences).
/// Identical to `hamming_distance`: both compare only sites where both sequences hold a base.
#[wasm_bindgen]
pub fn p_distance(seq1: &str, seq2: &str) -> f64 {
    hamming(&Mismatches::from_strs(seq1, seq2))
}

pub(crate) fn snps(m: &Mismatches) -> f64 {
    m.mismatches
}

pub(crate) fn hamming(m: &Mismatches) -> f64 {
    if m.sites == 0.0 {
        f64::NAN
    } else {
        m.mismatches / m.sites
    }
}
