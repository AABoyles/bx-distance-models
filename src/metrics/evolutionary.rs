/// Evolutionary distance models
use wasm_bindgen::prelude::*;
use crate::counts::{Counts, Mismatches};
use crate::metrics::basic::hamming;

/// Compute Tamura-Nei 93 (TN93) distance between two aligned sequences.
/// This is a more sophisticated nucleotide substitution model that accounts
/// for different rates of transitions and transversions.
#[wasm_bindgen]
pub fn tn93_distance(seq1: &str, seq2: &str) -> f64 {
    tn93(&Counts::from_strs(seq1, seq2))
}

/// Compute Jukes-Cantor (JC69) distance between two aligned sequences.
/// A simple evolutionary model assuming all substitutions are equally likely.
/// This is a basic but fundamental model used in phylogenetics.
#[wasm_bindgen]
pub fn jc69_distance(seq1: &str, seq2: &str) -> f64 {
    jc69(&Mismatches::from_strs(seq1, seq2))
}

/// Compute Kimura 2-Parameter (K2P) distance between two aligned sequences.
/// Accounts for different rates of transitions and transversions.
/// Simpler than TN93 but more accurate than Hamming or JC69.
#[wasm_bindgen]
pub fn k2p_distance(seq1: &str, seq2: &str) -> f64 {
    k2p(&Counts::from_strs(seq1, seq2))
}

/// Compute Tajima-Nei 92 (TN92) distance between two aligned sequences.
/// Similar to TN93 but with a slightly different parameterization.
#[wasm_bindgen]
pub fn tn92_distance(seq1: &str, seq2: &str) -> f64 {
    tn92(&Counts::from_strs(seq1, seq2))
}

pub(crate) fn tn93(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        return 1.0;
    }

    // Normalize frequencies
    let total2 = c.sites * 2.0;
    let g_a = c.freq[0] / total2;
    let g_c = c.freq[1] / total2;
    let g_g = c.freq[2] / total2;
    let g_t = c.freq[3] / total2;

    let g_r = g_a + g_g; // purine frequency
    let g_y = g_c + g_t; // pyrimidine frequency

    if g_r == 0.0 || g_y == 0.0 {
        return hamming(&c.into());
    }

    let p1 = c.ts_purine / c.sites;
    let p2 = c.ts_pyrimidine / c.sites;
    let q = c.transversions / c.sites;

    // TN93 formula
    let a1 = (2.0 * g_a * g_g) / g_r;
    let a2 = (2.0 * g_c * g_t) / g_y;
    let b = 2.0 * (g_r * g_y - (g_a * g_g * g_y) / g_r - (g_c * g_t * g_r) / g_y);

    // If A or G (C or T) never occurs, a1 (a2) and the matching transition proportion are both
    // 0, and that term of the formula vanishes: treat 0/0 and 0*ln as 0 rather than NaN.
    let ratio = |p: f64, a: f64| if a == 0.0 { 0.0 } else { p / a };
    let weighted_ln = |a: f64, term: f64| if a == 0.0 { 0.0 } else { a * term.ln() };

    let term1 = 1.0 - ratio(p1, a1) - q / (2.0 * g_r);
    let term2 = 1.0 - ratio(p2, a2) - q / (2.0 * g_y);
    let term3 = 1.0 - q / (2.0 * g_r * g_y);

    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 {
        // Distance is too large for the model; fall back to Hamming
        return hamming(&c.into());
    }

    let d = -weighted_ln(a1, term1) - weighted_ln(a2, term2) - b * term3.ln();

    d.max(0.0).min(1.0)
}

pub(crate) fn jc69(m: &Mismatches) -> f64 {
    if m.sites == 0.0 {
        return 1.0;
    }

    let p = m.mismatches / m.sites;

    // JC69 formula: d = -3/4 * ln(1 - 4/3 * p)
    let arg = 1.0 - (4.0 / 3.0) * p;
    if arg <= 0.0 {
        return 1.0;
    }

    (-0.75 * arg.ln()).min(1.0).max(0.0)
}

pub(crate) fn k2p(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        return 1.0;
    }

    let p = c.transitions() / c.sites;
    let q = c.transversions / c.sites;

    // K2P formula: d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
    let arg1 = 1.0 - 2.0 * p - q;
    let arg2 = 1.0 - 2.0 * q;

    if arg1 <= 0.0 || arg2 <= 0.0 {
        return 1.0;
    }

    let d = -0.5 * arg1.ln() - 0.25 * arg2.ln();
    d.min(1.0).max(0.0)
}

pub(crate) fn tn92(c: &Counts) -> f64 {
    if c.sites == 0.0 {
        return 1.0;
    }

    let total2 = c.sites * 2.0;
    let freq_a_f = c.freq[0] / total2;
    let freq_c_f = c.freq[1] / total2;
    let freq_g_f = c.freq[2] / total2;
    let freq_t_f = c.freq[3] / total2;

    let p = c.transitions() / c.sites;
    let q = c.transversions / c.sites;

    // TN92 formula
    let a_coeff = 2.0 * freq_a_f * freq_g_f;
    let c_coeff = 2.0 * freq_c_f * freq_t_f;
    let b_coeff = 2.0 * ((freq_a_f + freq_g_f) * (freq_c_f + freq_t_f));

    let term1 = 1.0 - p / a_coeff - q / (2.0 * (freq_a_f + freq_g_f));
    let term2 = 1.0 - p / c_coeff - q / (2.0 * (freq_c_f + freq_t_f));
    let term3 = 1.0 - q / b_coeff;

    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 {
        return hamming(&c.into());
    }

    let d = -a_coeff * term1.ln() - c_coeff * term2.ln() - (b_coeff / 2.0 - a_coeff - c_coeff) * term3.ln();
    d.max(0.0).min(1.0)
}
