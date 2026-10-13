//! Distances checked against R's `ape::dist.dna`, from the table that
//! `tests/reference/generate.R` writes (see that script for how the cases are built and where
//! ape's conventions differ from this library's).

use crate::counts::Counts;
use crate::encoding::Ambiguity;
use crate::metrics::{self, Metric};
use crate::sequence_set::SequenceSet;

const TABLE: &str = include_str!("../tests/reference/ape_distances.tsv");

/// ape model column → metrics that must reproduce it.
const MODELS: [(&str, &[Metric]); 8] = [
    ("N", &[Metric::Snps]),
    ("raw", &[Metric::Hamming, Metric::PDistance]),
    ("TS", &[Metric::Transitions]),
    ("TV", &[Metric::Transversions]),
    ("JC69", &[Metric::Jc69]),
    ("K80", &[Metric::K2p]),
    ("T92", &[Metric::Tn92]),
    ("TN93", &[Metric::Tn93]),
];

struct Case<'a> {
    id: &'a str,
    sites: f64,
    seq1: &'a str,
    seq2: &'a str,
    /// ape's value per column of `MODELS`; `None` where the table has NA (not comparable).
    expected: Vec<Option<f64>>,
}

fn cases() -> Vec<Case<'static>> {
    let mut lines = TABLE.lines().filter(|l| !l.starts_with('#'));
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let col = |name: &str| header.iter().position(|&h| h == name).unwrap_or_else(|| panic!("no column {name}"));
    let model_cols: Vec<usize> = MODELS.iter().map(|(name, _)| col(name)).collect();
    let (id, sites, seq1, seq2) = (col("id"), col("sites"), col("seq1"), col("seq2"));
    lines
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), header.len(), "malformed row: {line}");
            Case {
                id: f[id],
                sites: f[sites].parse().unwrap(),
                seq1: f[seq1],
                seq2: f[seq2],
                expected: model_cols
                    .iter()
                    .map(|&c| (f[c] != "NA").then(|| f[c].parse().unwrap_or_else(|_| panic!("bad number {}", f[c]))))
                    .collect(),
            }
        })
        .collect()
}

fn single_pair(metric: Metric, a: &str, b: &str) -> f64 {
    use metrics::{basic::*, diagnostics::*, evolutionary::*};
    match metric {
        Metric::Snps => snps_distance(a, b),
        Metric::Hamming => hamming_distance(a, b),
        Metric::PDistance => p_distance(a, b),
        Metric::Jc69 => jc69_distance(a, b),
        Metric::K2p => k2p_distance(a, b),
        Metric::Tn92 => tn92_distance(a, b),
        Metric::Tn93 => tn93_distance(a, b),
        Metric::Transitions => transitions_count(a, b),
        Metric::Transversions => transversions_count(a, b),
        Metric::TsTvRatio => ts_tv_ratio(a, b),
    }
}

/// Whether the library's `got` agrees with ape's `want` for a pair with `sites` compared sites.
fn agrees(got: f64, want: f64, sites: f64) -> bool {
    if want.is_finite() {
        // ape and this library evaluate the same formulas in a different order.
        (got - want).abs() <= 1e-12 * want.abs().max(1.0)
    } else if sites == 0.0 {
        got.is_nan()
    } else {
        // Saturated: ape takes the log of a non-positive number (NaN); the library gives Infinity.
        got == f64::INFINITY
    }
}

#[test]
fn reference_table_is_complete() {
    let cases = cases();
    assert!(cases.len() >= 200, "{} cases", cases.len());
    for (k, (name, _)) in MODELS.iter().enumerate() {
        let compared = cases.iter().filter(|c| c.expected[k].is_some()).count();
        let finite = cases.iter().filter(|c| c.expected[k].is_some_and(f64::is_finite)).count();
        assert!(compared >= 50 && finite >= 40, "{name}: {compared} compared, {finite} finite");
    }
    let saturated = cases.iter().filter(|c| c.sites > 0.0 && c.expected[5].is_some_and(f64::is_nan)).count();
    assert!(saturated >= 5, "only {saturated} pairs saturate K80");
}

#[test]
fn compared_sites_match_ape() {
    for c in cases() {
        assert_eq!(Counts::from_strs(c.seq1, c.seq2).sites, c.sites, "case {}", c.id);
    }
}

#[test]
fn string_functions_match_ape() {
    for c in cases() {
        for (k, (name, metrics)) in MODELS.iter().enumerate() {
            let Some(want) = c.expected[k] else { continue };
            for &metric in *metrics {
                let got = single_pair(metric, c.seq1, c.seq2);
                assert!(agrees(got, want, c.sites), "case {} {name} {metric:?}: got {got}, ape {want}", c.id);
            }
        }
    }
}

#[test]
fn sequence_set_matches_ape() {
    let cases = cases();
    let seqs: Vec<&str> = cases.iter().flat_map(|c| [c.seq1, c.seq2]).collect();
    let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
    for (n, c) in cases.iter().enumerate() {
        for (k, (name, metrics)) in MODELS.iter().enumerate() {
            let Some(want) = c.expected[k] else { continue };
            for &metric in *metrics {
                let got = set.pair(metric, 2 * n, 2 * n + 1);
                assert!(agrees(got, want, c.sites), "case {} {name} {metric:?}: got {got}, ape {want}", c.id);
            }
        }
    }
}
