mod encoding;
mod counts;
mod metrics;
mod pairwise;
mod sequence_set;

#[cfg(test)]
mod tests;

use wasm_bindgen::prelude::*;

#[cfg(feature = "console_error_panic_hook")]
use console_error_panic_hook;

/// Initialize panic hook for better error messages in browser console
#[wasm_bindgen(start)]
pub fn init() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

// Re-export all public functions for WASM bindings
pub use metrics::basic::*;
pub use metrics::evolutionary::*;
pub use metrics::diagnostics::*;
pub use pairwise::*;
pub use sequence_set::{Links, SequenceSet};
