/// Module declarations and re-exports for metrics
pub mod basic;
pub mod evolutionary;
pub mod diagnostics;

#[allow(unused_imports)]
pub use basic::{snps_distance, hamming_distance, p_distance};
#[allow(unused_imports)]
pub use evolutionary::{tn93_distance, jc69_distance, k2p_distance, tn92_distance};
#[allow(unused_imports)]
pub use diagnostics::{transitions_count, transversions_count, ts_tv_ratio};

use crate::counts::{Counts, CountsWithGc, Mismatches};

/// What a per-pair scan collects, cheapest first. A batch over several metrics uses the largest
/// scan any of them needs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Scan {
    /// Compared sites and mismatches ([`Mismatches`], 2 popcounts per 64 sites).
    Mismatches,
    /// Full [`Counts`] (~10 popcounts per 64 sites), without per-sequence G+C content.
    Counts,
    /// Full counts plus the first sequence's G+C content ([`CountsWithGc`]); T92 only.
    CountsWithGc,
}

/// Every metric, evaluated from the per-pair [`Counts`], or from [`Mismatches`] when that suffices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Metric {
    Snps,
    Hamming,
    PDistance,
    Jc69,
    K2p,
    Tn92,
    Tn93,
    Transitions,
    Transversions,
    TsTvRatio,
}

pub const METRIC_NAMES: &str =
    "snps, hamming, p_distance, jc69, k2p, tn92, tn93, transitions, transversions, ts_tv_ratio";

impl Metric {
    pub fn parse(name: &str) -> Option<Metric> {
        Some(match name.to_ascii_lowercase().as_str() {
            "snps" | "snp" => Metric::Snps,
            "hamming" => Metric::Hamming,
            "p" | "p_distance" | "pdistance" => Metric::PDistance,
            "jc69" => Metric::Jc69,
            "k2p" | "k80" => Metric::K2p,
            "tn92" => Metric::Tn92,
            "tn93" => Metric::Tn93,
            "transitions" => Metric::Transitions,
            "transversions" => Metric::Transversions,
            "ts_tv_ratio" | "tstv" | "ts_tv" => Metric::TsTvRatio,
            _ => return None,
        })
    }

    /// The cheapest per-pair scan that provides everything the metric needs.
    pub fn scan(self) -> Scan {
        match self {
            Metric::Snps | Metric::Hamming | Metric::PDistance | Metric::Jc69 => Scan::Mismatches,
            Metric::Tn92 => Scan::CountsWithGc,
            _ => Scan::Counts,
        }
    }

    // `eval` and `eval_mismatches` are kept out of line. `f64::ln` is an LLVM intrinsic that may
    // be speculated, so when the match was inlined into the dense row loop LLVM hoisted the
    // formulas' `ln` calls above it and every pair paid for every metric's logarithms
    // (measured: `transitions` +90 ns/pair, Hamming +25% of its time).
    #[inline(never)]
    pub fn eval(self, c: &Counts) -> f64 {
        match self {
            Metric::Snps | Metric::Hamming | Metric::PDistance | Metric::Jc69 => self.eval_mismatches(&c.into()),
            Metric::K2p => evolutionary::k2p(c),
            Metric::Tn92 => unreachable!("T92 needs per-sequence G+C content: use eval_with_gc"),
            Metric::Tn93 => evolutionary::tn93(c),
            Metric::Transitions => diagnostics::transitions(c),
            Metric::Transversions => diagnostics::transversions(c),
            Metric::TsTvRatio => diagnostics::ts_tv(c),
        }
    }

    /// Evaluates any metric, including T92 (the only one that needs `gc1`).
    #[inline(never)]
    pub fn eval_with_gc(self, c: &CountsWithGc) -> f64 {
        match self {
            Metric::Tn92 => evolutionary::tn92(c),
            _ => self.eval(&c.counts),
        }
    }

    /// Evaluates a metric whose [`scan`](Self::scan) is [`Scan::Mismatches`].
    #[inline(never)]
    pub fn eval_mismatches(self, m: &Mismatches) -> f64 {
        match self {
            Metric::Snps => basic::snps(m),
            Metric::Hamming | Metric::PDistance => basic::hamming(m),
            Metric::Jc69 => evolutionary::jc69(m),
            _ => unreachable!("{self:?} needs full counts"),
        }
    }

    /// Mismatch count above which a pair can no longer score `<= threshold`, given that it can
    /// compare at most `max_sites` sites. Used to stop scanning a pair early; exact, not a
    /// heuristic: the SNP count is the mismatch count, and Hamming, JC69, K2P, T92 and TN93 are
    /// never smaller than the p-distance (mismatches / sites; applying -ln(1 - x) >= x to each log
    /// term of the corrected models gives exactly p). Other metrics are never pruned.
    pub fn mismatch_bound(self, threshold: f64, max_sites: u32) -> f64 {
        match self {
            Metric::Snps => threshold,
            Metric::Hamming
            | Metric::PDistance
            | Metric::Jc69
            | Metric::K2p
            | Metric::Tn92
            | Metric::Tn93 => threshold * max_sites as f64,
            _ => f64::INFINITY,
        }
    }
}
