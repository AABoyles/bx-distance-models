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

use crate::counts::Counts;

/// Every metric, evaluated from the same per-pair [`Counts`].
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

    pub fn eval(self, c: &Counts) -> f64 {
        match self {
            Metric::Snps => basic::snps(c),
            Metric::Hamming | Metric::PDistance => basic::hamming(c),
            Metric::Jc69 => evolutionary::jc69(c),
            Metric::K2p => evolutionary::k2p(c),
            Metric::Tn92 => evolutionary::tn92(c),
            Metric::Tn93 => evolutionary::tn93(c),
            Metric::Transitions => diagnostics::transitions(c),
            Metric::Transversions => diagnostics::transversions(c),
            Metric::TsTvRatio => diagnostics::ts_tv(c),
        }
    }

    /// Mismatch count above which a pair can no longer score `<= threshold`, given that it can
    /// compare at most `max_sites` sites. Used to stop scanning a pair early; exact, not a
    /// heuristic: the SNP count is the mismatch count, and Hamming, JC69, K2P and TN93 are never
    /// smaller than the p-distance (mismatches / sites). Other metrics are never pruned.
    pub fn mismatch_bound(self, threshold: f64, max_sites: u32) -> f64 {
        match self {
            Metric::Snps => threshold,
            Metric::Hamming | Metric::PDistance | Metric::Jc69 | Metric::K2p | Metric::Tn93 => {
                threshold * max_sites as f64
            }
            _ => f64::INFINITY,
        }
    }
}
