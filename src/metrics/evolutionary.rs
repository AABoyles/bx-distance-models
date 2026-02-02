/// Evolutionary distance models
use wasm_bindgen::prelude::*;
use crate::utils::*;
use crate::metrics::basic::hamming_distance;

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

/// Compute Jukes-Cantor (JC69) distance between two aligned sequences.
/// A simple evolutionary model assuming all substitutions are equally likely.
/// This is a basic but fundamental model used in phylogenetics.
#[wasm_bindgen]
pub fn jc69_distance(seq1: &str, seq2: &str) -> f64 {
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
        return 1.0;
    }

    let p = mismatches as f64 / compared as f64;

    // JC69 formula: d = -3/4 * ln(1 - 4/3 * p)
    let arg = 1.0 - (4.0 / 3.0) * p;
    if arg <= 0.0 {
        return 1.0;
    }

    (-0.75 * arg.ln()).min(1.0).max(0.0)
}

/// Compute Kimura 2-Parameter (K2P) distance between two aligned sequences.
/// Accounts for different rates of transitions and transversions.
/// Simpler than TN93 but more accurate than Hamming or JC69.
#[wasm_bindgen]
pub fn k2p_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 1.0;
    }

    let mut transitions = 0u32;
    let mut transversions = 0u32;
    let mut compared = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        // Skip gaps and N's
        if a == b'-' || b == b'-' || a == b'N' || b == b'N' {
            continue;
        }

        if a == b {
            compared += 1;
            continue;
        }

        compared += 1;

        // Classify substitution type
        if is_purine(a) && is_purine(b) {
            transitions += 1;
        } else if is_pyrimidine(a) && is_pyrimidine(b) {
            transitions += 1;
        } else {
            transversions += 1;
        }
    }

    if compared == 0 {
        return 1.0;
    }

    let p = transitions as f64 / compared as f64;
    let q = transversions as f64 / compared as f64;

    // K2P formula: d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
    let arg1 = 1.0 - 2.0 * p - q;
    let arg2 = 1.0 - 2.0 * q;

    if arg1 <= 0.0 || arg2 <= 0.0 {
        return 1.0;
    }

    let d = -0.5 * arg1.ln() - 0.25 * arg2.ln();
    d.min(1.0).max(0.0)
}

/// Compute Tajima-Nei 92 (TN92) distance between two aligned sequences.
/// Similar to TN93 but with a slightly different parameterization.
#[wasm_bindgen]
pub fn tn92_distance(seq1: &str, seq2: &str) -> f64 {
    let bytes1 = seq1.as_bytes();
    let bytes2 = seq2.as_bytes();
    let len = bytes1.len().min(bytes2.len());

    if len == 0 {
        return 1.0;
    }

    let mut total_sites = 0u32;
    let mut transitions = 0u32;
    let mut transversions = 0u32;
    let mut freq_a = 0u32;
    let mut freq_c = 0u32;
    let mut freq_g = 0u32;
    let mut freq_t = 0u32;

    for i in 0..len {
        let a = bytes1[i];
        let b = bytes2[i];

        if !is_valid_base(a) || !is_valid_base(b) {
            continue;
        }

        total_sites += 1;

        // Count frequencies from both sequences
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

        if is_purine(a) && is_purine(b) {
            transitions += 1;
        } else if is_pyrimidine(a) && is_pyrimidine(b) {
            transitions += 1;
        } else {
            transversions += 1;
        }
    }

    if total_sites == 0 {
        return 1.0;
    }

    let total2 = (total_sites * 2) as f64;
    let freq_a_f = freq_a as f64 / total2;
    let freq_c_f = freq_c as f64 / total2;
    let freq_g_f = freq_g as f64 / total2;
    let freq_t_f = freq_t as f64 / total2;

    let p = transitions as f64 / total_sites as f64;
    let q = transversions as f64 / total_sites as f64;

    // TN92 formula
    let a_coeff = 2.0 * freq_a_f * freq_g_f;
    let c_coeff = 2.0 * freq_c_f * freq_t_f;
    let b_coeff = 2.0 * ((freq_a_f + freq_g_f) * (freq_c_f + freq_t_f));

    let term1 = 1.0 - p / a_coeff - q / (2.0 * (freq_a_f + freq_g_f));
    let term2 = 1.0 - p / c_coeff - q / (2.0 * (freq_c_f + freq_t_f));
    let term3 = 1.0 - q / b_coeff;

    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 {
        return hamming_distance(seq1, seq2);
    }

    let d = -a_coeff * term1.ln() - c_coeff * term2.ln() - (b_coeff / 2.0 - a_coeff - c_coeff) * term3.ln();
    d.max(0.0).min(1.0)
}
