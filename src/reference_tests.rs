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

// ---------------------------------------------------------------------------------------------
// HIV-TRACE's tn93, from the table that `tests/reference/tn93/generate.py` writes: every pair of
// the FASTA files in `tests/reference/tn93/`, in every ambiguity mode.

mod tn93 {
    use std::collections::HashMap;

    use crate::counts::Counts;
    use crate::encoding::Ambiguity;
    use crate::metrics::{evolutionary::tn93_distance, Metric};
    use crate::sequence_set::SequenceSet;

    const TABLE: &str = include_str!("../tests/reference/tn93_distances.tsv");

    const FILES: [(&str, &str); 5] = [
        ("distances.fasta", include_str!("../tests/reference/tn93/distances.fasta")),
        ("ambiguity.fasta", include_str!("../tests/reference/tn93/ambiguity.fasta")),
        ("hivtrace-g.fasta", include_str!("../tests/reference/tn93/hivtrace-g.fasta")),
        ("n_codes.fasta", include_str!("../tests/reference/tn93/n_codes.fasta")),
        ("random.fasta", include_str!("../tests/reference/tn93/random.fasta")),
    ];

    /// `(name, sequence)` in file order; sequences may span several lines.
    fn fasta(text: &str) -> Vec<(&str, String)> {
        let mut out: Vec<(&str, String)> = Vec::new();
        for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
            match line.strip_prefix('>') {
                Some(name) => out.push((name, String::new())),
                None => out.last_mut().expect("sequence before a name").1.push_str(line),
            }
        }
        out
    }

    struct Row<'a> {
        file: &'a str,
        mode: &'a str,
        fraction: Option<f64>,
        id1: &'a str,
        id2: &'a str,
        /// `None`: tn93 left the pair out (at most one compared site).
        tn93: Option<f64>,
    }

    fn rows() -> Vec<Row<'static>> {
        let mut lines = TABLE.lines().filter(|l| !l.starts_with('#'));
        assert_eq!(lines.next(), Some("file\tmode\tfraction\tid1\tid2\ttn93"));
        lines
            .map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                assert_eq!(f.len(), 6, "malformed row: {line}");
                Row {
                    file: f[0],
                    mode: f[1],
                    fraction: (f[2] != "-").then(|| f[2].parse().unwrap()),
                    id1: f[3],
                    id2: f[4],
                    tn93: (f[5] != "NA").then(|| f[5].parse().unwrap()),
                }
            })
            .collect()
    }

    /// A FASTA file encoded for one configuration.
    struct Loaded {
        set: SequenceSet,
        index: HashMap<&'static str, usize>,
        seqs: Vec<String>,
    }

    /// One [`Loaded`] set per file and configuration, built on first use.
    #[derive(Default)]
    struct Sets(HashMap<(&'static str, &'static str), Loaded>);

    impl Sets {
        fn get(&mut self, file: &'static str, mode: &'static str) -> &Loaded {
            self.0.entry((file, mode)).or_insert_with(|| {
                let text = FILES.iter().find(|(name, _)| *name == file).unwrap_or_else(|| panic!("no file {file}")).1;
                let records = fasta(text);
                let index = records.iter().enumerate().map(|(k, (name, _))| (*name, k)).collect();
                let seqs: Vec<String> = records.into_iter().map(|(_, s)| s).collect();
                let mode = Ambiguity::parse(mode).unwrap_or_else(|| panic!("unknown mode {mode}"));
                Loaded { set: SequenceSet::from_strs(&seqs, mode), index, seqs }
            })
        }
    }

    /// Whether the library's `got` agrees with tn93's `want`.
    fn agrees(got: f64, want: f64) -> bool {
        if want == 1000.0 {
            // tn93's TN93_MAX_DIST: the model saturated
            got == f64::INFINITY
        } else {
            (got - want).abs() <= 1e-12 * want.abs().max(1.0)
        }
    }

    /// Whether the library implements tn93's handling of this configuration yet.
    fn supported(fraction: Option<f64>) -> bool {
        fraction.is_none()
    }

    #[test]
    fn table_covers_every_mode() {
        let rows = rows();
        for mode in ["skip", "resolve", "average", "gapmm"] {
            let finite = rows.iter().filter(|r| r.mode == mode && r.tn93.is_some_and(|d| d < 1000.0)).count();
            assert!(finite >= 300, "{mode}: {finite} finite values");
        }
    }

    #[test]
    fn sequence_set_matches_tn93() {
        let mut sets = Sets::default();
        let (mut checked, mut missing_base) = (0, 0);
        for row in rows() {
            let Loaded { set, index, seqs } = sets.get(row.file, row.mode);
            let (i, j) = (index[row.id1], index[row.id2]);
            if !supported(row.fraction) {
                continue;
            }
            let id = format!("{} -a {} {} {}", row.file, row.mode, row.id1, row.id2);
            let counts = set.tally::<Counts>(i, j, f64::INFINITY).unwrap();
            let Some(want) = row.tn93 else {
                assert!(counts.sites <= 1.0, "{id}: tn93 left out a pair with {} sites", counts.sites);
                continue;
            };
            if counts.freq.contains(&0.0) {
                // tn93 falls back to K2P when a base is missing from the pair
                missing_base += 1;
                continue;
            }
            let got = set.pair(Metric::Tn93, i, j);
            assert!(agrees(got, want), "{id}: got {got}, tn93 {want}");
            if row.mode == "skip" {
                let single = tn93_distance(&seqs[i], &seqs[j]);
                assert!(agrees(single, want), "{id}: tn93_distance gives {single}, tn93 {want}");
            }
            checked += 1;
        }
        assert!(checked >= 1400 && missing_base <= 10, "{checked} pairs checked, {missing_base} missing a base");
    }

    /// `links()` stops scanning a pair early once it cannot pass the threshold; on tn93's own
    /// values, it must keep exactly the pairs within it.
    #[test]
    fn links_keep_the_pairs_tn93_keeps() {
        let mut sets = Sets::default();
        let rows = rows();
        let mut groups: Vec<(&str, &str, Option<f64>)> = rows.iter().map(|r| (r.file, r.mode, r.fraction)).collect();
        groups.dedup();
        for (file, mode, fraction) in groups {
            if !supported(fraction) {
                continue;
            }
            let Loaded { set, index, .. } = sets.get(file, mode);
            for threshold in [0.0, 0.01, 0.03, 0.1, 0.5] {
                let mut want: Vec<(u32, u32)> = rows
                    .iter()
                    .filter(|r| (r.file, r.mode, r.fraction) == (file, mode, fraction))
                    .map(|r| (index[r.id1], index[r.id2], r.tn93))
                    .filter(|&(_, _, d)| d.is_some_and(|d| d <= threshold))
                    .map(|(i, j, _)| (i as u32, j as u32))
                    .collect();
                want.sort();
                let links = set.sparse(Metric::Tn93, threshold, 0..set.len());
                let got: Vec<(u32, u32)> = links.source.iter().zip(&links.target).map(|(&i, &j)| (i, j)).collect();
                assert_eq!(got, want, "{file} -a {mode} -g {fraction:?}, threshold {threshold}");
            }
        }
    }
}
