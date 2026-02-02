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
