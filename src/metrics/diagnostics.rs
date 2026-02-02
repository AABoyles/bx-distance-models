/// Transition/transversion diagnostics
use wasm_bindgen::prelude::*;
use crate::utils::*;

/// Count the number of transitions between two aligned sequences.
/// Transitions are purine-to-purine (A<->G) or pyrimidine-to-pyrimidine (C<->T) changes.
#[wasm_bindgen]
pub fn transitions_count(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    let mut count = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        if a == b'-' || b == b'-' || a == b'N' || b == b'N' {
            continue;
        }

        if a == b {
            continue;
        }

        if (is_purine(a) && is_purine(b)) || (is_pyrimidine(a) && is_pyrimidine(b)) {
            count += 1;
        }
    }

    count as f64
}

/// Count the number of transversions between two aligned sequences.
/// Transversions are purine-to-pyrimidine or pyrimidine-to-purine changes.
#[wasm_bindgen]
pub fn transversions_count(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    let mut count = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        if a == b'-' || b == b'-' || a == b'N' || b == b'N' {
            continue;
        }

        if a == b {
            continue;
        }

        if (is_purine(a) && is_pyrimidine(b)) || (is_pyrimidine(a) && is_purine(b)) {
            count += 1;
        }
    }

    count as f64
}

/// Compute the transition/transversion ratio (Ts/Tv) for two sequences.
/// Returns the ratio of transitions to transversions. Returns 0 if no transversions.
#[wasm_bindgen]
pub fn ts_tv_ratio(seq1: &str, seq2: &str) -> f64 {
    let transitions = transitions_count(seq1, seq2) as i32;
    let transversions = transversions_count(seq1, seq2) as i32;

    if transversions == 0 {
        if transitions == 0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        transitions as f64 / transversions as f64
    }
}
