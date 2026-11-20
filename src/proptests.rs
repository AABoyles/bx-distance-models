//! Property-based tests: the batch kernels against a per-site reference, and invariants every
//! metric must satisfy, on generated sets of related sequences (failures shrink to a minimal set).
//!
//! `PROPTEST_CASES=10000 cargo test proptests` runs a longer search than the default 256 cases.

use proptest::prelude::*;

use crate::counts::{Counts, CountsWithGc, Mismatches};
use crate::encoding::{Ambiguity, IUPAC};
use crate::metrics::{self, Metric};
use crate::sequence_set::{row_splits, Links, SequenceSet, SetOptions};

const ALL: [Metric; 10] = [
    Metric::Snps, Metric::Hamming, Metric::PDistance, Metric::Jc69, Metric::K2p,
    Metric::Tn92, Metric::Tn93, Metric::Transitions, Metric::Transversions, Metric::TsTvRatio,
];
const MODES: [Ambiguity; 4] = [Ambiguity::Skip, Ambiguity::Resolve, Ambiguity::Average, Ambiguity::GapMm];
/// Characters substituted into sequences: bases in either case, U, gaps, N, every partial
/// ambiguity code and characters that are not nucleotides at all.
const NOISE: &[u8] = b"ACGTacgtUu--NnRYSWKMBDHVrykm?.*X";

/// A root sequence of 0-300 bases, so lengths straddle the 64-site words.
fn root() -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(prop::sample::select(b"ACGT".to_vec()), 0..300)
}

/// 2-7 sequences derived from one root: substitutions, noise, and sometimes a shorter length.
fn sequence_set() -> impl Strategy<Value = Vec<String>> {
    let derived = (
        prop::collection::vec((any::<prop::sample::Index>(), prop::sample::select(b"ACGT".to_vec())), 0..40),
        prop::collection::vec((any::<prop::sample::Index>(), prop::sample::select(NOISE.to_vec())), 0..20),
        prop::option::weighted(0.2, any::<prop::sample::Index>()),
    );
    (root(), prop::collection::vec(derived, 2..8)).prop_map(|(root, edits)| {
        edits
            .into_iter()
            .map(|(subs, noise, cut)| {
                let mut s = root.clone();
                if !s.is_empty() {
                    for (at, c) in subs.into_iter().chain(noise) {
                        s[at.index(root.len())] = c;
                    }
                }
                if let Some(cut) = cut {
                    s.truncate(cut.index(root.len() + 1));
                }
                String::from_utf8(s).unwrap()
            })
            .collect()
    })
}

/// Per-site reference: one IUPAC lookup per site, no bit planes, no sparse ambiguity lists.
/// Follows tn93's `computeTN93`: N is an ambiguity code where the mode compares it, and in GapMM
/// mode a gap facing anything but a gap is an N.
fn reference_counts(a: &str, b: &str, mode: Ambiguity) -> CountsWithGc {
    let mut c = CountsWithGc::default();
    for (&x, &y) in a.as_bytes().iter().zip(b.as_bytes()) {
        let (mut cx, mut cy) = (IUPAC[x as usize], IUPAC[y as usize]);
        if !mode.compares_n() {
            cx = if cx == 15 { 0 } else { cx };
            cy = if cy == 15 { 0 } else { cy };
        }
        if mode.gap_is_n() && (cx != 0 || cy != 0) {
            cx = if cx == 0 { 15 } else { cx };
            cy = if cy == 0 { 15 } else { cy };
        }
        if cx == 0 || cy == 0 {
            continue;
        }
        if cx.count_ones() == 1 && cy.count_ones() == 1 {
            c.add(cx.trailing_zeros() as usize, cy.trailing_zeros() as usize, 1.0);
        } else if mode != Ambiguity::Skip {
            c.add_codes(cx, cy, 1.0, mode);
        }
    }
    c
}

/// tn93 `-a resolve -g fraction` for one pair, per site: resolved, unless the sites resolution
/// touched (`ambig_count` in `computeTN93`: a base facing a code that contains it, or two
/// ambiguity codes) are at least `fraction` of the compared sites; then averaged.
fn reference_max_ambiguity(a: &str, b: &str, fraction: f64) -> CountsWithGc {
    let (mut sites, mut ambiguous) = (0u32, 0u32);
    for (&x, &y) in a.as_bytes().iter().zip(b.as_bytes()) {
        let (cx, cy) = (IUPAC[x as usize], IUPAC[y as usize]);
        if cx == 0 || cy == 0 {
            continue;
        }
        sites += 1;
        let (ax, ay) = (cx.count_ones() > 1, cy.count_ones() > 1);
        ambiguous += ((ax || ay) && (cx & cy != 0 || (ax && ay))) as u32;
    }
    let too_ambiguous = fraction < 1.0 && sites as f64 * fraction <= ambiguous as f64;
    reference_counts(a, b, if too_ambiguous { Ambiguity::Average } else { Ambiguity::Resolve })
}

/// Every configuration a set can be built with: each mode, and resolve with a maximum
/// ambiguity fraction.
fn config() -> impl Strategy<Value = (Ambiguity, Option<f64>)> {
    prop_oneof![
        prop::sample::select(MODES.to_vec()).prop_map(|m| (m, None)),
        prop_oneof![Just(0.0), Just(0.015), 0.0..0.2f64, Just(1.0)].prop_map(|f| (Ambiguity::Resolve, Some(f))),
    ]
}

fn near(x: f64, y: f64) -> bool {
    (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits() || (x - y).abs() <= 1e-9 * (1.0 + x.abs().max(y.abs()))
}

fn near_counts(x: &CountsWithGc, y: &CountsWithGc) -> bool {
    let (a, b) = (&x.counts, &y.counts);
    near(a.sites, b.sites) && near(a.ts_purine, b.ts_purine) && near(a.ts_pyrimidine, b.ts_pyrimidine)
        && near(a.transversions, b.transversions) && (0..4).all(|k| near(a.freq[k], b.freq[k]))
        && near(x.gc1, y.gc1)
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

fn is_distance(metric: Metric) -> bool {
    matches!(metric, Metric::Hamming | Metric::PDistance | Metric::Jc69 | Metric::K2p | Metric::Tn92 | Metric::Tn93)
}

proptest! {
    #![proptest_config(ProptestConfig { failure_persistence: None, ..ProptestConfig::default() })]

    /// Every scan level, in every ambiguity mode, counts what a per-site loop counts; in skip
    /// mode, so do the string kernels behind the single-pair functions, bit for bit.
    #[test]
    fn kernels_match_per_site_reference(seqs in sequence_set()) {
        for mode in MODES {
            let set = SequenceSet::from_strs(&seqs, mode);
            for i in 0..seqs.len() {
                for j in 0..seqs.len() {
                    let (a, b) = (&seqs[i], &seqs[j]);
                    let want = reference_counts(a, b, mode);
                    let full = set.tally::<CountsWithGc>(i, j, f64::INFINITY).unwrap();
                    prop_assert!(near_counts(&full, &want), "{mode:?} {i},{j}\n{a}\n{b}\n{full:?}\n{want:?}");
                    prop_assert_eq!(set.tally::<Counts>(i, j, f64::INFINITY).unwrap(), full.counts);
                    let light = set.tally::<Mismatches>(i, j, f64::INFINITY).unwrap();
                    prop_assert!(near(light.sites, full.counts.sites) && near(light.mismatches, full.counts.mismatches()));
                    if mode == Ambiguity::Skip {
                        prop_assert_eq!(full, want);
                        prop_assert_eq!(CountsWithGc::from_strs(a, b), want);
                        prop_assert_eq!(Mismatches::from_strs(a, b), Mismatches::from(&want.counts));
                    }
                }
            }
        }
    }

    /// The batch API and the single-pair functions give the same value for every metric.
    #[test]
    fn set_and_single_pair_functions_agree(seqs in sequence_set()) {
        let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
        for i in 0..seqs.len() {
            for j in 0..seqs.len() {
                for metric in ALL {
                    let (got, want) = (set.pair(metric, i, j), single_pair(metric, &seqs[i], &seqs[j]));
                    prop_assert!(got.to_bits() == want.to_bits(), "{metric:?} {i},{j}: {got} vs {want}");
                }
            }
        }
    }

    /// The counts are symmetric, and with integer counts so is every metric, bit for bit. With
    /// ambiguity resolution the counts are fractional and summed in a different order, so a model
    /// near saturation (a log argument near zero) can amplify their last-bit differences.
    #[test]
    fn metrics_are_symmetric(seqs in sequence_set()) {
        for mode in MODES {
            let set = SequenceSet::from_strs(&seqs, mode);
            for i in 0..seqs.len() {
                for j in i + 1..seqs.len() {
                    let ij = set.tally::<CountsWithGc>(i, j, f64::INFINITY).unwrap();
                    let ji = set.tally::<CountsWithGc>(j, i, f64::INFINITY).unwrap();
                    let gc = ij.counts.freq[1] + ij.counts.freq[2];
                    let swapped = CountsWithGc { counts: ji.counts, gc1: gc - ji.gc1 };
                    prop_assert!(near_counts(&ij, &swapped), "{mode:?} {i},{j}\n{ij:?}\n{ji:?}");
                    if mode == Ambiguity::Skip {
                        for metric in ALL {
                            let (d_ij, d_ji) = (set.pair(metric, i, j), set.pair(metric, j, i));
                            prop_assert!(d_ij.to_bits() == d_ji.to_bits(), "{metric:?} {i},{j}: {d_ij} vs {d_ji}");
                        }
                    }
                }
            }
        }
    }

    /// A sequence is at distance 0 from itself, or NaN (undefined) without a comparable site.
    /// Average mode is the exception for ambiguity codes: R vs R averages A/G mismatches in.
    #[test]
    fn a_sequence_is_at_distance_zero_from_itself(seqs in sequence_set()) {
        for mode in [Ambiguity::Skip, Ambiguity::Resolve] {
            let set = SequenceSet::from_strs(&seqs, mode);
            for (i, seq) in seqs.iter().enumerate() {
                let sites = set.tally::<Counts>(i, i, f64::INFINITY).unwrap().sites;
                for metric in ALL {
                    let d = set.pair(metric, i, i);
                    let ok = if is_distance(metric) && sites == 0.0 { d.is_nan() } else { d == 0.0 };
                    prop_assert!(ok, "{mode:?} {metric:?} seq {i} ({sites} sites): {d}\n{seq}");
                }
            }
        }
    }

    /// Ranges, bounds and relations that hold for every pair.
    #[test]
    fn values_are_in_range(seqs in sequence_set()) {
        for mode in MODES {
            let set = SequenceSet::from_strs(&seqs, mode);
            for i in 0..seqs.len() {
                for j in i + 1..seqs.len() {
                    let c = set.tally::<Counts>(i, j, f64::INFINITY).unwrap();
                    let d = |m| set.pair(m, i, j);
                    let p = d(Metric::PDistance);
                    if c.sites == 0.0 {
                        for metric in ALL.into_iter().filter(|&m| is_distance(m)) {
                            prop_assert!(d(metric).is_nan(), "{metric:?}");
                        }
                        continue;
                    }
                    prop_assert!((0.0..=1.0).contains(&p));
                    prop_assert!(d(Metric::Snps) <= c.sites + 1e-9);
                    prop_assert!(near(d(Metric::Transitions) + d(Metric::Transversions), d(Metric::Snps)));
                    // The corrected distances are never below the p-distance (links() prunes on it).
                    for metric in [Metric::Jc69, Metric::K2p, Metric::Tn92, Metric::Tn93] {
                        prop_assert!(d(metric) >= p, "{mode:?} {metric:?} {} < p {p}\n{}\n{}", d(metric), seqs[i], seqs[j]);
                    }
                }
            }
        }
    }

    /// Dense output in any number of row ranges, any combination of metrics, joins up to the
    /// full triangle, and each value is the pair's own. A batch scans each pair once at the level
    /// its most demanding metric needs, so with ambiguity resolution (fractional counts) a light
    /// metric can differ from its own scan in the last bits; with integer counts it is exact.
    #[test]
    fn dense_rows_join_up_to_the_triangle(
        seqs in sequence_set(),
        picks in prop::collection::vec(prop::sample::select(ALL.to_vec()), 1..4),
        parts in 1usize..6,
        mode in prop::sample::select(MODES.to_vec()),
    ) {
        let n = seqs.len();
        let set = SequenceSet::from_strs(&seqs, mode);
        let full = set.dense(&picks, 0..n);
        let np = n * (n - 1) / 2;
        let mut k = 0;
        for i in 0..n {
            for j in i + 1..n {
                for (m, &metric) in picks.iter().enumerate() {
                    let (got, want) = (full[m * np + k], set.pair(metric, i, j));
                    let exact = got.to_bits() == want.to_bits();
                    prop_assert!(exact || (mode != Ambiguity::Skip && near(got, want)), "{mode:?} {metric:?} {i},{j}: {got} vs {want}");
                }
                k += 1;
            }
        }
        let splits = row_splits(n, parts);
        let mut joined = vec![Vec::new(); picks.len()];
        for w in splits.windows(2) {
            let rows = w[0] as usize..w[1] as usize;
            let part = set.dense(&picks, rows);
            let width = part.len() / picks.len();
            for (m, out) in joined.iter_mut().enumerate() {
                out.extend_from_slice(&part[m * width..(m + 1) * width]);
            }
        }
        let joined: Vec<f64> = joined.concat();
        prop_assert_eq!(joined.len(), full.len());
        prop_assert!(joined.iter().zip(&full).all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    /// links() is exactly the dense output filtered by the threshold, early exit and all.
    #[test]
    fn links_are_filtered_dense_output(
        seqs in sequence_set(),
        metric in prop::sample::select(ALL.to_vec()),
        threshold in prop_oneof![0.0..0.3f64, 0.0..20.0f64, Just(0.0), Just(f64::INFINITY)],
        (mode, fraction) in config(),
        start in any::<prop::sample::Index>(),
    ) {
        let n = seqs.len();
        let options = SetOptions { max_ambiguity_fraction: fraction };
        let set = SequenceSet::with_options(&seqs, mode, options).unwrap();
        let rows = start.index(n + 1)..n;
        let dense = set.dense(&[metric], rows.clone());
        let mut want = Links::default();
        let mut k = 0;
        for i in rows.clone() {
            for j in i + 1..n {
                if dense[k] <= threshold {
                    want.source.push(i as u32);
                    want.target.push(j as u32);
                    want.distance.push(dense[k]);
                }
                k += 1;
            }
        }
        prop_assert_eq!(set.sparse(metric, threshold, rows), want);
    }

    /// With a maximum ambiguity fraction, each pair is what the per-site reference resolves or
    /// averages, as tn93 decides, at every scan level.
    #[test]
    fn max_ambiguity_fraction_resolves_or_averages_like_tn93(
        seqs in sequence_set(),
        fraction in prop_oneof![Just(0.0), Just(0.015), 0.0..0.2f64, Just(1.0)],
    ) {
        let options = SetOptions { max_ambiguity_fraction: Some(fraction) };
        let set = SequenceSet::with_options(&seqs, Ambiguity::Resolve, options).unwrap();
        for i in 0..seqs.len() {
            for j in 0..seqs.len() {
                let (a, b) = (&seqs[i], &seqs[j]);
                let want = reference_max_ambiguity(a, b, fraction);
                let full = set.tally::<CountsWithGc>(i, j, f64::INFINITY).unwrap();
                prop_assert!(near_counts(&full, &want), "-g {fraction} {i},{j}\n{a}\n{b}\n{full:?}\n{want:?}");
                prop_assert_eq!(set.tally::<Counts>(i, j, f64::INFINITY).unwrap(), full.counts);
                let light = set.tally::<Mismatches>(i, j, f64::INFINITY).unwrap();
                prop_assert!(near(light.sites, full.counts.sites) && near(light.mismatches, full.counts.mismatches()));
            }
        }
    }
}
