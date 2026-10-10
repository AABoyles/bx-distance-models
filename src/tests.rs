#[cfg(test)]
mod tests {
    use crate::counts::Counts;
    use crate::encoding::{Ambiguity, IUPAC};
    use crate::metrics::{self, Metric};
    use crate::sequence_set::{pairs_in_rows, row_splits, Links, SequenceSet};

    #[test]
    fn test_hamming_identical() {
        assert_eq!(metrics::basic::hamming_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_hamming_different() {
        assert_eq!(metrics::basic::hamming_distance("AAAA", "TTTT"), 1.0);
    }

    #[test]
    fn test_hamming_partial() {
        assert_eq!(metrics::basic::hamming_distance("AAGT", "ACGT"), 0.25);
    }

    #[test]
    fn test_hamming_with_gaps() {
        // With gap, positions: A vs A (same), - vs C (skip), G vs G (same), T vs T (same)
        // Only 3 positions compared, 0 mismatches
        assert_eq!(metrics::basic::hamming_distance("A-GT", "ACGT"), 0.0);
    }

    #[test]
    fn test_hamming_with_mismatch_and_gap() {
        // With gap, positions: A vs A (same), - vs C (skip), T vs T (same), T vs T (same)
        // 3 positions compared, 0 mismatches
        assert_eq!(metrics::basic::hamming_distance("A-TT", "ACTT"), 0.0);
    }

    #[test]
    fn test_tn93_identical() {
        assert_eq!(metrics::evolutionary::tn93_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_snps_identical() {
        assert_eq!(metrics::basic::snps_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_snps_different() {
        assert_eq!(metrics::basic::snps_distance("AAAA", "TTTT"), 4.0);
    }

    #[test]
    fn test_snps_with_gaps() {
        // With gap, positions: A vs A (skip), - vs C (skip), G vs T (diff), T vs T (same)
        // 2 positions compared, 1 mismatch
        assert_eq!(metrics::basic::snps_distance("A-GT", "ACTT"), 1.0);
    }

    #[test]
    fn test_jc69_identical() {
        assert_eq!(metrics::evolutionary::jc69_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_k2p_identical() {
        assert_eq!(metrics::evolutionary::k2p_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_tn92_identical() {
        assert_eq!(metrics::evolutionary::tn92_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_p_distance_identical() {
        assert_eq!(metrics::basic::p_distance("ACGT", "ACGT"), 0.0);
    }

    #[test]
    fn test_transitions_count() {
        // A->G is a transition
        assert_eq!(metrics::diagnostics::transitions_count("AA", "GG"), 2.0);
    }

    #[test]
    fn test_transversions_count() {
        // A->T is a transversion
        assert_eq!(metrics::diagnostics::transversions_count("AA", "TT"), 2.0);
    }

    #[test]
    fn test_ts_tv_ratio() {
        // 2 transitions, 2 transversions = ratio of 1.0
        let seq1 = "AAAT";
        let seq2 = "GGTG";
        let ratio = metrics::diagnostics::ts_tv_ratio(seq1, seq2);
        assert_eq!(ratio, 1.0);
    }

    // ------------------------------------------------------------------
    // Equivalence with the previous byte-loop implementation
    // ------------------------------------------------------------------

    /// Verbatim copies of the implementations this crate shipped before the encoded kernel,
    /// used as the reference for inputs where both agree on which sites count (ACGTN-).
    mod legacy {
        fn is_valid_base(b: u8) -> bool { matches!(b, b'A' | b'C' | b'G' | b'T') }
        fn is_purine(b: u8) -> bool { matches!(b, b'A' | b'G') }
        fn is_pyrimidine(b: u8) -> bool { matches!(b, b'C' | b'T') }
        fn skip(a: u8, b: u8) -> bool { a == b'-' || b == b'-' || a == b'N' || b == b'N' }

        pub fn snps(s1: &str, s2: &str) -> f64 {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            if len == 0 { return 0.0; }
            let mut mismatches = 0u32;
            for i in 0..len {
                if skip(b1[i], b2[i]) { continue; }
                if b1[i] != b2[i] { mismatches += 1; }
            }
            mismatches as f64
        }

        pub fn hamming(s1: &str, s2: &str) -> f64 {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            if len == 0 { return 1.0; }
            let (mut mismatches, mut compared) = (0u32, 0u32);
            for i in 0..len {
                if skip(b1[i], b2[i]) { continue; }
                compared += 1;
                if b1[i] != b2[i] { mismatches += 1; }
            }
            if compared == 0 { 1.0 } else { mismatches as f64 / compared as f64 }
        }

        pub fn p_distance(s1: &str, s2: &str) -> f64 {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            if len == 0 { return 1.0; }
            let (mut differences, mut compared) = (0u32, 0u32);
            for i in 0..len {
                if !is_valid_base(b1[i]) || !is_valid_base(b2[i]) { continue; }
                compared += 1;
                if b1[i] != b2[i] { differences += 1; }
            }
            if compared == 0 { return 1.0; }
            differences as f64 / compared as f64
        }

        fn class_counts(s1: &str, s2: &str) -> (u32, u32, u32, u32, [u32; 4]) {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            let (mut total, mut ts1, mut ts2, mut tv) = (0u32, 0u32, 0u32, 0u32);
            let mut f = [0u32; 4];
            for i in 0..len {
                let (a, b) = (b1[i], b2[i]);
                if !is_valid_base(a) || !is_valid_base(b) { continue; }
                total += 1;
                for x in [a, b] {
                    match x { b'A' => f[0] += 1, b'C' => f[1] += 1, b'G' => f[2] += 1, b'T' => f[3] += 1, _ => {} }
                }
                if a == b { continue; }
                if is_purine(a) && is_purine(b) { ts1 += 1; }
                else if is_pyrimidine(a) && is_pyrimidine(b) { ts2 += 1; }
                else { tv += 1; }
            }
            (total, ts1, ts2, tv, f)
        }

        pub fn tn93(s1: &str, s2: &str) -> f64 {
            if s1.is_empty() || s2.is_empty() { return 1.0; }
            let (total_sites, transitions1, transitions2, transversions, f) = class_counts(s1, s2);
            if total_sites == 0 { return 1.0; }
            let total2 = (total_sites * 2) as f64;
            let (g_a, g_c, g_g, g_t) = (f[0] as f64 / total2, f[1] as f64 / total2, f[2] as f64 / total2, f[3] as f64 / total2);
            let g_r = g_a + g_g;
            let g_y = g_c + g_t;
            if g_r == 0.0 || g_y == 0.0 { return hamming(s1, s2); }
            let p1 = transitions1 as f64 / total_sites as f64;
            let p2 = transitions2 as f64 / total_sites as f64;
            let q = transversions as f64 / total_sites as f64;
            let a1 = (2.0 * g_a * g_g) / g_r;
            let a2 = (2.0 * g_c * g_t) / g_y;
            let b = 2.0 * (g_r * g_y - (g_a * g_g * g_y) / g_r - (g_c * g_t * g_r) / g_y);
            let term1 = 1.0 - p1 / a1 - q / (2.0 * g_r);
            let term2 = 1.0 - p2 / a2 - q / (2.0 * g_y);
            let term3 = 1.0 - q / (2.0 * g_r * g_y);
            if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 { return hamming(s1, s2); }
            let d = -a1 * term1.ln() - a2 * term2.ln() - b * term3.ln();
            d.max(0.0).min(1.0)
        }

        pub fn tn92(s1: &str, s2: &str) -> f64 {
            if s1.is_empty() || s2.is_empty() { return 1.0; }
            let (total_sites, ts1, ts2, transversions, f) = class_counts(s1, s2);
            if total_sites == 0 { return 1.0; }
            let total2 = (total_sites * 2) as f64;
            let (fa, fc, fg, ft) = (f[0] as f64 / total2, f[1] as f64 / total2, f[2] as f64 / total2, f[3] as f64 / total2);
            let p = (ts1 + ts2) as f64 / total_sites as f64;
            let q = transversions as f64 / total_sites as f64;
            let a_coeff = 2.0 * fa * fg;
            let c_coeff = 2.0 * fc * ft;
            let b_coeff = 2.0 * ((fa + fg) * (fc + ft));
            let term1 = 1.0 - p / a_coeff - q / (2.0 * (fa + fg));
            let term2 = 1.0 - p / c_coeff - q / (2.0 * (fc + ft));
            let term3 = 1.0 - q / b_coeff;
            if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 { return hamming(s1, s2); }
            let d = -a_coeff * term1.ln() - c_coeff * term2.ln() - (b_coeff / 2.0 - a_coeff - c_coeff) * term3.ln();
            d.max(0.0).min(1.0)
        }

        pub fn jc69(s1: &str, s2: &str) -> f64 {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            if len == 0 { return 1.0; }
            let (mut mismatches, mut compared) = (0u32, 0u32);
            for i in 0..len {
                if skip(b1[i], b2[i]) { continue; }
                compared += 1;
                if b1[i] != b2[i] { mismatches += 1; }
            }
            if compared == 0 { return 1.0; }
            let p = mismatches as f64 / compared as f64;
            let arg = 1.0 - (4.0 / 3.0) * p;
            if arg <= 0.0 { return 1.0; }
            (-0.75 * arg.ln()).min(1.0).max(0.0)
        }

        /// (transitions, transversions, compared) with the '-'/'N' skip rule of K2P and the
        /// diagnostics.
        fn ts_tv_counts(s1: &str, s2: &str) -> (u32, u32, u32) {
            let (b1, b2) = (s1.as_bytes(), s2.as_bytes());
            let len = b1.len().min(b2.len());
            let (mut ts, mut tv, mut compared) = (0u32, 0u32, 0u32);
            for i in 0..len {
                let (a, b) = (b1[i], b2[i]);
                if skip(a, b) { continue; }
                compared += 1;
                if a == b { continue; }
                if (is_purine(a) && is_purine(b)) || (is_pyrimidine(a) && is_pyrimidine(b)) { ts += 1; }
                else { tv += 1; }
            }
            (ts, tv, compared)
        }

        pub fn k2p(s1: &str, s2: &str) -> f64 {
            if s1.is_empty() || s2.is_empty() { return 1.0; }
            let (transitions, transversions, compared) = ts_tv_counts(s1, s2);
            if compared == 0 { return 1.0; }
            let p = transitions as f64 / compared as f64;
            let q = transversions as f64 / compared as f64;
            let arg1 = 1.0 - 2.0 * p - q;
            let arg2 = 1.0 - 2.0 * q;
            if arg1 <= 0.0 || arg2 <= 0.0 { return 1.0; }
            let d = -0.5 * arg1.ln() - 0.25 * arg2.ln();
            d.min(1.0).max(0.0)
        }

        pub fn transitions(s1: &str, s2: &str) -> f64 {
            ts_tv_counts(s1, s2).0 as f64
        }

        pub fn transversions(s1: &str, s2: &str) -> f64 {
            ts_tv_counts(s1, s2).1 as f64
        }

        pub fn ts_tv(s1: &str, s2: &str) -> f64 {
            let (ts, tv) = (transitions(s1, s2), transversions(s1, s2));
            if tv == 0.0 { if ts == 0.0 { 0.0 } else { f64::INFINITY } } else { ts / tv }
        }
    }

    /// Deterministic xorshift generator so failures reproduce.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
        fn chance(&mut self, p: f64) -> bool { (self.next() % 1_000_000) as f64 / 1_000_000.0 < p }
    }

    /// Pairs of related sequences over `alphabet`, lengths straddling 64-site word boundaries,
    /// divergence from identical to saturated, sometimes of unequal length.
    fn random_pairs(seed: u64, count: usize, alphabet: &[u8]) -> Vec<(String, String)> {
        let mut rng = Rng(seed);
        let rates = [0.0, 0.01, 0.05, 0.2, 0.5, 0.9];
        (0..count)
            .map(|k| {
                let len = 40 + rng.below(300);
                let a: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
                let rate = rates[k % rates.len()];
                let mut b: Vec<u8> = a.iter().map(|&c| if rng.chance(rate) { b"ACGT"[rng.below(4)] } else { c }).collect();
                let mut a = a;
                for seq in [&mut a, &mut b] {
                    for c in seq.iter_mut() {
                        if rng.chance(0.05) { *c = alphabet[rng.below(alphabet.len())]; }
                    }
                }
                if k % 5 == 0 { b.truncate(len - rng.below(len / 2)); }
                (String::from_utf8(a).unwrap(), String::from_utf8(b).unwrap())
            })
            .collect()
    }

    fn same(a: f64, b: f64) -> bool {
        a.to_bits() == b.to_bits()
    }

    const ALL: [Metric; 10] = [
        Metric::Snps, Metric::Hamming, Metric::PDistance, Metric::Jc69, Metric::K2p,
        Metric::Tn92, Metric::Tn93, Metric::Transitions, Metric::Transversions, Metric::TsTvRatio,
    ];

    fn legacy_eval(metric: Metric, a: &str, b: &str) -> f64 {
        match metric {
            Metric::Snps => legacy::snps(a, b),
            Metric::Hamming => legacy::hamming(a, b),
            Metric::PDistance => legacy::p_distance(a, b),
            Metric::Jc69 => legacy::jc69(a, b),
            Metric::K2p => legacy::k2p(a, b),
            Metric::Tn92 => legacy::tn92(a, b),
            Metric::Tn93 => legacy::tn93(a, b),
            Metric::Transitions => legacy::transitions(a, b),
            Metric::Transversions => legacy::transversions(a, b),
            Metric::TsTvRatio => legacy::ts_tv(a, b),
        }
    }

    #[test]
    fn matches_legacy_implementation_bit_for_bit() {
        let pairs = random_pairs(0x9E3779B97F4A7C15, 600, b"ACGTN-");
        for (a, b) in &pairs {
            let set = SequenceSet::from_strs(&[a, b], Ambiguity::Skip);
            let counts = Counts::from_strs(a, b);
            for metric in ALL {
                let want = legacy_eval(metric, a, b);
                assert!(same(metric.eval(&counts), want), "{metric:?} strings\n{a}\n{b}");
                assert!(same(set.pair(metric, 0, 1), want), "{metric:?} SequenceSet\n{a}\n{b}");
            }
        }
    }

    #[test]
    fn string_functions_match_legacy() {
        for (a, b) in random_pairs(7, 100, b"ACGTN-") {
            assert!(same(metrics::basic::snps_distance(&a, &b), legacy::snps(&a, &b)));
            assert!(same(metrics::basic::hamming_distance(&a, &b), legacy::hamming(&a, &b)));
            assert!(same(metrics::basic::p_distance(&a, &b), legacy::p_distance(&a, &b)));
            assert!(same(metrics::evolutionary::tn93_distance(&a, &b), legacy::tn93(&a, &b)));
            assert!(same(metrics::evolutionary::tn92_distance(&a, &b), legacy::tn92(&a, &b)));
            assert!(same(metrics::evolutionary::jc69_distance(&a, &b), legacy::jc69(&a, &b)));
            assert!(same(metrics::evolutionary::k2p_distance(&a, &b), legacy::k2p(&a, &b)));
            assert!(same(metrics::diagnostics::ts_tv_ratio(&a, &b), legacy::ts_tv(&a, &b)));
        }
    }

    // ------------------------------------------------------------------
    // Character handling
    // ------------------------------------------------------------------

    #[test]
    fn lower_case_and_u_are_bases() {
        assert_eq!(metrics::basic::hamming_distance("acgu", "ACGT"), 0.0);
        assert_eq!(metrics::evolutionary::tn93_distance("acgtacgt", "ACGTACGT"), 0.0);
        assert_eq!(metrics::basic::snps_distance("acgt", "ACGA"), 1.0);
    }

    #[test]
    fn ambiguity_codes_and_other_characters_are_skipped_by_default() {
        // R, ?, . and N are not compared: 3 sites, no mismatches
        assert_eq!(metrics::basic::hamming_distance("ARGT?.N", "AAGTACG"), 0.0);
        assert_eq!(Counts::from_strs("ARGT?.N", "AAGTACG").sites, 3.0);
    }

    #[test]
    fn tn93_handles_an_absent_base() {
        // No G anywhere: the purine-transition term vanishes instead of turning into NaN
        // (which the old implementation clamped to 0.0).
        let (a, b) = ("AACCTTAACCTT", "AACCTTAACCTA");
        assert_eq!(legacy::tn93(a, b), 0.0);
        let d = metrics::evolutionary::tn93_distance(a, b);
        assert!(d.is_finite() && d >= metrics::basic::hamming_distance(a, b), "{d}");
    }

    // ------------------------------------------------------------------
    // Ambiguity fallback
    // ------------------------------------------------------------------

    fn set_counts(a: &str, b: &str, mode: Ambiguity) -> Counts {
        SequenceSet::from_strs(&[a, b], mode).counts(0, 1, f64::INFINITY).unwrap()
    }

    #[test]
    fn resolve_mode_counts_compatible_codes_as_matches() {
        let c = set_counts("ARGT", "AAGT", Ambiguity::Resolve);
        assert_eq!((c.sites, c.mismatches()), (4.0, 0.0));
        assert_eq!(c.freq, [4.0, 0.0, 2.0, 2.0]);
        // R vs R shares A and G: a match split between them
        let c = set_counts("R", "R", Ambiguity::Resolve);
        assert_eq!((c.sites, c.mismatches(), c.freq), (1.0, 0.0, [1.0, 0.0, 1.0, 0.0]));
    }

    #[test]
    fn resolve_mode_averages_incompatible_codes() {
        // R = A|G vs C: A-C and G-C, both transversions
        let c = set_counts("ARGT", "ACGT", Ambiguity::Resolve);
        assert_eq!((c.sites, c.transversions, c.transitions()), (4.0, 1.0, 0.0));
    }

    #[test]
    fn average_mode_spreads_every_resolution() {
        // R vs A: half A-A (match), half G-A (transition)
        let c = set_counts("R", "A", Ambiguity::Average);
        assert_eq!((c.sites, c.ts_purine, c.freq), (1.0, 0.5, [1.5, 0.0, 0.5, 0.0]));
        // R vs R: AA, AG, GA, GG
        let c = set_counts("R", "R", Ambiguity::Average);
        assert_eq!((c.sites, c.ts_purine), (1.0, 0.5));
        // Y vs R: every resolution is a transversion
        let c = set_counts("Y", "R", Ambiguity::Average);
        assert_eq!((c.sites, c.transversions), (1.0, 1.0));
        // B (C|G|T) vs A: three transversions-or-transitions, weighted 1/3 each
        let c = set_counts("B", "A", Ambiguity::Average);
        assert!((c.ts_purine - 1.0 / 3.0).abs() < 1e-12 && (c.transversions - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn n_and_gaps_are_never_compared() {
        for mode in [Ambiguity::Skip, Ambiguity::Resolve, Ambiguity::Average] {
            let c = set_counts("ANR-T", "ARNRT", mode);
            assert_eq!(c.sites, 2.0, "{mode:?}");
        }
    }

    #[test]
    fn skip_mode_ignores_ambiguous_sites() {
        for (a, b) in random_pairs(11, 100, b"ACGTN-RYKMSWBDHV") {
            let set = SequenceSet::from_strs(&[&a, &b], Ambiguity::Skip);
            assert_eq!(set.ambiguous_count(), 0);
            assert_eq!(set.counts(0, 1, f64::INFINITY).unwrap(), Counts::from_strs(&a, &b));
        }
    }

    /// Straightforward per-site reference for the fast kernel + ambiguity fallback.
    fn reference_counts(a: &str, b: &str, mode: Ambiguity) -> Counts {
        let mut c = Counts::default();
        for (&x, &y) in a.as_bytes().iter().zip(b.as_bytes()) {
            let (cx, cy) = (IUPAC[x as usize], IUPAC[y as usize]);
            let informative = |code: u8| code != 0 && code != 15;
            if !informative(cx) || !informative(cy) { continue; }
            let single = |code: u8| code.count_ones() == 1;
            if single(cx) && single(cy) {
                c.add(cx.trailing_zeros() as usize, cy.trailing_zeros() as usize, 1.0);
            } else if mode != Ambiguity::Skip {
                c.add_codes(cx, cy, 1.0, mode);
            }
        }
        c
    }

    fn close(a: &Counts, b: &Counts) -> bool {
        let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * (1.0 + x.abs());
        near(a.sites, b.sites) && near(a.ts_purine, b.ts_purine) && near(a.ts_pyrimidine, b.ts_pyrimidine)
            && near(a.transversions, b.transversions) && (0..4).all(|k| near(a.freq[k], b.freq[k]))
    }

    #[test]
    fn ambiguity_fallback_matches_per_site_reference() {
        for mode in [Ambiguity::Resolve, Ambiguity::Average] {
            for (a, b) in random_pairs(23, 300, b"ACGTN-RYKMSWBDHVacgtry") {
                let got = set_counts(&a, &b, mode);
                let want = reference_counts(&a, &b, mode);
                assert!(close(&got, &want), "{mode:?}\n{a}\n{b}\n{got:?}\n{want:?}");
            }
        }
    }

    #[test]
    fn unambiguous_pairs_use_the_fast_path_in_every_mode() {
        let seqs = ["ACGTACGTAC", "ACGTTCGTAC", "ACRTACGTAC"];
        let set = SequenceSet::from_strs(&seqs, Ambiguity::Resolve);
        assert_eq!(set.ambiguous_count(), 1);
        assert_eq!(set.counts(0, 1, f64::INFINITY).unwrap(), Counts::from_strs(seqs[0], seqs[1]));
    }

    // ------------------------------------------------------------------
    // Batch output: dense, sparse, row ranges
    // ------------------------------------------------------------------

    fn random_set(seed: u64, n: usize, alphabet: &[u8], mode: Ambiguity) -> (Vec<String>, SequenceSet) {
        let mut rng = Rng(seed);
        let len = 150 + rng.below(100);
        let root: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
        let seqs: Vec<String> = (0..n)
            .map(|k| {
                let rate = [0.002, 0.01, 0.03, 0.1, 0.4][k % 5];
                root.iter()
                    .map(|&c| if rng.chance(rate) { b"ACGT"[rng.below(4)] } else if rng.chance(0.02) { alphabet[rng.below(alphabet.len())] } else { c })
                    .map(char::from)
                    .collect()
            })
            .collect();
        let set = SequenceSet::from_strs(&seqs, mode);
        (seqs, set)
    }

    #[test]
    fn dense_matches_pairwise_in_triangle_order() {
        let (_, set) = random_set(3, 25, b"ACGTN-", Ambiguity::Skip);
        let dense = set.dense(&[Metric::Tn93, Metric::Snps], 0..set.len());
        let np = 25 * 24 / 2;
        assert_eq!(dense.len(), 2 * np);
        for i in 0..25 {
            for j in i + 1..25 {
                let k = crate::pairwise::get_distance_index(i as u32, j as u32, 25) as usize;
                assert!(same(dense[k], set.pair(Metric::Tn93, i, j)));
                assert!(same(dense[np + k], set.pair(Metric::Snps, i, j)));
            }
        }
    }

    #[test]
    fn row_ranges_concatenate_to_the_full_triangle() {
        let (_, set) = random_set(5, 37, b"ACGTN-", Ambiguity::Skip);
        let full = set.dense(&[Metric::Hamming], 0..37);
        for parts in [1, 2, 3, 4, 7, 50] {
            let splits = row_splits(37, parts);
            assert_eq!(splits.len(), parts + 1);
            assert_eq!((splits[0], splits[parts]), (0, 37));
            assert!(splits.windows(2).all(|w| w[0] <= w[1]));
            let joined: Vec<f64> = splits.windows(2)
                .flat_map(|w| set.dense(&[Metric::Hamming], w[0] as usize..w[1] as usize))
                .collect();
            assert_eq!(joined, full, "parts={parts}");
        }
    }

    #[test]
    fn row_splits_balance_pairs() {
        let splits = row_splits(1000, 4);
        let sizes: Vec<u64> = splits.windows(2).map(|w| pairs_in_rows(1000, w[0] as usize..w[1] as usize)).collect();
        let (min, max) = (*sizes.iter().min().unwrap(), *sizes.iter().max().unwrap());
        assert!(max - min < 2000, "{sizes:?}");
    }

    #[test]
    fn links_equal_filtered_dense_output_with_pruning() {
        let cases = [
            (Metric::Tn93, 0.02), (Metric::Tn93, 0.15), (Metric::Hamming, 0.05), (Metric::Jc69, 0.05),
            (Metric::K2p, 0.05), (Metric::Snps, 3.0), (Metric::Tn92, 0.05), (Metric::Transitions, 2.0),
        ];
        let mut total_links = 0;
        for mode in [Ambiguity::Skip, Ambiguity::Resolve, Ambiguity::Average] {
            let (_, set) = random_set(17, 30, b"ACGTN-RY", mode);
            for (metric, threshold) in cases {
                for rows in [0..30, 4..19] {
                    let dense = set.dense(&[metric], rows.clone());
                    let mut want = Links::default();
                    let mut k = 0;
                    for i in rows.clone() {
                        for j in i + 1..30 {
                            if dense[k] <= threshold {
                                want.source.push(i as u32);
                                want.target.push(j as u32);
                                want.distance.push(dense[k]);
                            }
                            k += 1;
                        }
                    }
                    let got = set.sparse(metric, threshold, rows.clone());
                    assert_eq!(got, want, "{mode:?} {metric:?} {threshold} {rows:?}");
                    total_links += want.distance.len();
                }
            }
        }
        assert!(total_links > 500, "test data should produce links: {total_links}");
    }

    #[test]
    fn pruning_bound_is_exact_on_random_pairs() {
        // TN93 / JC69 / K2P are never below the p-distance, which is what pruning relies on.
        for (a, b) in random_pairs(31, 600, b"ACGTN-") {
            let c = Counts::from_strs(&a, &b);
            let p = metrics::basic::hamming(&c);
            for metric in [Metric::Tn93, Metric::Jc69, Metric::K2p] {
                assert!(metric.eval(&c) >= p, "{metric:?}\n{a}\n{b}");
            }
        }
    }

    #[test]
    fn small_and_empty_sets() {
        for n in [0, 1] {
            let seqs = vec!["ACGT"; n];
            let set = SequenceSet::from_strs(&seqs, Ambiguity::Skip);
            assert!(set.dense(&[Metric::Tn93], 0..n).is_empty());
            assert_eq!(set.sparse(Metric::Tn93, 1.0, 0..n), Links::default());
            let splits = row_splits(n, 3);
            assert_eq!((splits.len(), splits[0], splits[3]), (4, 0, n as u32));
        }
    }

    #[test]
    fn distance_index_has_no_intermediate_overflow() {
        // i * (2n - i - 1) exceeds u32 here, the index itself does not
        let n = 80_000u32;
        let (i, j) = (n - 2, n - 1);
        let expected = (n as u64 * (n as u64 - 1) / 2 - 1) as u32;
        assert_eq!(crate::pairwise::get_distance_index(i, j, n), expected);
    }

    #[test]
    fn metric_names_parse() {
        assert_eq!(Metric::parse("TN93"), Some(Metric::Tn93));
        assert_eq!(Metric::parse("p_distance"), Some(Metric::PDistance));
        assert_eq!(Metric::parse("ts_tv_ratio"), Some(Metric::TsTvRatio));
        assert_eq!(Metric::parse("nope"), None);
        assert_eq!(Ambiguity::parse("Average"), Some(Ambiguity::Average));
        assert_eq!(Ambiguity::parse("ignore"), None);
    }
}
