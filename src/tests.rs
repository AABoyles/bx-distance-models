#[cfg(test)]
mod tests {
    use crate::metrics;

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
}
