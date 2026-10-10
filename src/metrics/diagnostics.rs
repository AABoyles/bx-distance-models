/// Transition/transversion diagnostics
use wasm_bindgen::prelude::*;
use crate::counts::Counts;

/// Count the number of transitions between two aligned sequences.
/// Transitions are purine-to-purine (A<->G) or pyrimidine-to-pyrimidine (C<->T) changes.
#[wasm_bindgen]
pub fn transitions_count(seq1: &str, seq2: &str) -> f64 {
    transitions(&Counts::from_strs(seq1, seq2))
}

/// Count the number of transversions between two aligned sequences.
/// Transversions are purine-to-pyrimidine or pyrimidine-to-purine changes.
#[wasm_bindgen]
pub fn transversions_count(seq1: &str, seq2: &str) -> f64 {
    transversions(&Counts::from_strs(seq1, seq2))
}

/// Compute the transition/transversion ratio (Ts/Tv) for two sequences.
/// Returns the ratio of transitions to transversions; 0 if there are neither,
/// infinity if there are transitions but no transversions.
#[wasm_bindgen]
pub fn ts_tv_ratio(seq1: &str, seq2: &str) -> f64 {
    ts_tv(&Counts::from_strs(seq1, seq2))
}

pub(crate) fn transitions(c: &Counts) -> f64 {
    c.transitions()
}

pub(crate) fn transversions(c: &Counts) -> f64 {
    c.transversions
}

pub(crate) fn ts_tv(c: &Counts) -> f64 {
    let (ts, tv) = (c.transitions(), c.transversions);
    if tv == 0.0 {
        if ts == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        ts / tv
    }
}
