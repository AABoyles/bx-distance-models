/// Evolutionary distance models
///
/// Each returns the expected number of substitutions per site under its model. These distances
/// are not bounded by 1. When the sequences are too divergent for a model (one of its logarithms
/// would need a non-positive argument), the distance is `f64::INFINITY`.
use wasm_bindgen::prelude::*;
use crate::counts::{Counts, CountsWithGc, Mismatches};

/// Compute Tamura-Nei 93 (TN93) distance between two aligned sequences.
/// This is a more sophisticated nucleotide substitution model that accounts
/// for different rates of transitions and transversions.
/// Substitutions per site (may exceed 1); `Infinity` when the model saturates.
#[wasm_bindgen]
pub fn tn93_distance(seq1: &str, seq2: &str) -> f64 {
    tn93(&Counts::from_strs(seq1, seq2))
}

/// Compute Jukes-Cantor (JC69) distance between two aligned sequences.
/// A simple evolutionary model assuming all substitutions are equally likely.
/// This is a basic but fundamental model used in phylogenetics.
/// Substitutions per site (may exceed 1); `Infinity` when p >= 3/4.
#[wasm_bindgen]
pub fn jc69_distance(seq1: &str, seq2: &str) -> f64 {
    jc69(&Mismatches::from_strs(seq1, seq2))
}

/// Compute Kimura 2-Parameter (K2P) distance between two aligned sequences.
/// Accounts for different rates of transitions and transversions.
/// Simpler than TN93 but more accurate than Hamming or JC69.
/// Substitutions per site (may exceed 1); `Infinity` when the model saturates.
#[wasm_bindgen]
pub fn k2p_distance(seq1: &str, seq2: &str) -> f64 {
    k2p(&Counts::from_strs(seq1, seq2))
}

/// Compute the Tamura (1992) distance (T92) between two aligned sequences.
/// Kimura's two-parameter model extended for G+C-content bias, using the G+C content of each
/// sequence. Exposed as `tn92` for compatibility.
/// Substitutions per site (may exceed 1); `Infinity` when the model saturates.
#[wasm_bindgen]
pub fn tn92_distance(seq1: &str, seq2: &str) -> f64 {
    tn92(&CountsWithGc::from_strs(seq1, seq2))
}

/// `num / den`, taking a 0/0 term of a degenerate base composition as 0, its limit.
fn ratio(num: f64, den: f64) -> f64 {
    if den == 0.0 { 0.0 } else { num / den }
}

/// `weight * ln(term)`, taken as 0 when the weight is 0 (the term drops out of the model).
fn weighted_ln(weight: f64, term: f64) -> f64 {
    if weight == 0.0 { 0.0 } else { weight * term.ln() }
}

/// Whether a log term that is part of the model needs a non-positive argument.
fn saturated(weight: f64, term: f64) -> bool {
    weight != 0.0 && term <= 0.0
}

pub(crate) fn tn93(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        return f64::NAN;
    }

    // Normalize frequencies
    let total2 = c.sites * 2.0;
    let g_a = c.freq[0] / total2;
    let g_c = c.freq[1] / total2;
    let g_g = c.freq[2] / total2;
    let g_t = c.freq[3] / total2;

    let g_r = g_a + g_g; // purine frequency
    let g_y = g_c + g_t; // pyrimidine frequency

    let p1 = c.ts_purine / c.sites;
    let p2 = c.ts_pyrimidine / c.sites;
    let q = c.transversions / c.sites;

    // TN93 formula. If a base, or every purine or every pyrimidine, never occurs, the terms that
    // involve it are 0/0 (no substitution of that kind is possible either); they drop out.
    let a1 = ratio(2.0 * g_a * g_g, g_r);
    let a2 = ratio(2.0 * g_c * g_t, g_y);
    let b = 2.0 * (g_r * g_y - ratio(g_a * g_g * g_y, g_r) - ratio(g_c * g_t * g_r, g_y));

    let term1 = 1.0 - ratio(p1, a1) - ratio(q, 2.0 * g_r);
    let term2 = 1.0 - ratio(p2, a2) - ratio(q, 2.0 * g_y);
    let term3 = 1.0 - ratio(q, 2.0 * g_r * g_y);

    if saturated(a1, term1) || saturated(a2, term2) || saturated(b, term3) {
        return f64::INFINITY;
    }

    let d = -weighted_ln(a1, term1) - weighted_ln(a2, term2) - weighted_ln(b, term3);
    d.max(0.0)
}

pub(crate) fn jc69(m: &Mismatches) -> f64 {
    if m.sites == 0.0 {
        return f64::NAN;
    }

    let p = m.mismatches / m.sites;

    // JC69 formula: d = -3/4 * ln(1 - 4/3 * p)
    let arg = 1.0 - (4.0 / 3.0) * p;
    if arg <= 0.0 {
        return f64::INFINITY;
    }

    (-0.75 * arg.ln()).max(0.0)
}

pub(crate) fn k2p(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        return f64::NAN;
    }

    let p = c.transitions() / c.sites;
    let q = c.transversions / c.sites;

    // K2P formula: d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
    let arg1 = 1.0 - 2.0 * p - q;
    let arg2 = 1.0 - 2.0 * q;

    if arg1 <= 0.0 || arg2 <= 0.0 {
        return f64::INFINITY;
    }

    let d = -0.5 * arg1.ln() - 0.25 * arg2.ln();
    d.max(0.0)
}

pub(crate) fn tn92(gc: &CountsWithGc) -> f64 {
    let c = &gc.counts;
    if c.sites == 0.0 {
        return f64::NAN;
    }

    let p = c.transitions() / c.sites;
    let q = c.transversions / c.sites;

    // G+C content of each sequence at the compared sites
    let theta1 = gc.gc1 / c.sites;
    let theta2 = (c.freq[1] + c.freq[2] - gc.gc1) / c.sites;

    // Tamura (1992): d = -h ln(1 - P/h - Q) - (1 - h)/2 ln(1 - 2Q), h = θ1 + θ2 - 2 θ1 θ2.
    // h = 0 means both sequences are all A/T (or all G/C), so no transition is possible and the
    // first term drops out.
    let h = theta1 + theta2 - 2.0 * theta1 * theta2;
    let term1 = 1.0 - ratio(p, h) - q;
    let term2 = 1.0 - 2.0 * q;

    if saturated(h, term1) || saturated(1.0 - h, term2) {
        return f64::INFINITY;
    }

    let d = -weighted_ln(h, term1) - weighted_ln(0.5 * (1.0 - h), term2);
    d.max(0.0)
}
