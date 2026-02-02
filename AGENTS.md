# AGENTS.md - Project Architecture Guide

## Quick Facts
- **Type**: Rust library compiled to WebAssembly (WASM)
- **Language**: Rust (src/) + JavaScript/TypeScript (bindings)
- **Build Tool**: wasm-pack → generates JS/WASM in `pkg/`
- **Target**: Browser & Node.js via WASM bindings
- **Purpose**: Fast nucleotide sequence distance metrics

## Architecture Overview

```
bx-distance-models/
├── src/
│   ├── lib.rs                    # Main entry + re-exports for WASM
│   ├── utils.rs                  # Base validation, nucleotide classification
│   ├── pairwise.rs               # Batch distance computation (JS interop)
│   ├── metrics/
│   │   ├── mod.rs                # Module declarations
│   │   ├── basic.rs              # SNPs, Hamming, P-distance
│   │   ├── evolutionary.rs       # TN93, JC69, K2P, TN92
│   │   └── diagnostics.rs        # Transitions/transversions/ratios
│   └── tests.rs                  # Unit tests
├── pkg/                          # Generated WASM bindings (build output)
├── Cargo.toml                    # Rust dependencies + build config
├── package.json                  # NPM metadata, build scripts
└── target/                       # Rust build artifacts
```

## Module Responsibilities

### `lib.rs` - Main Entry Point
- Declares all submodules
- Re-exports public functions for WASM binding via `pub use`
- Calls `console_error_panic_hook::set_once()` on init for better browser debugging
- All public functions are automatically available in JavaScript

### `utils.rs` - Utility Functions
- `is_valid_base(byte)` → checks if A/T/C/G/N
- `is_purine(byte)` → A or G
- `is_pyrimidine(byte)` → C or T
- Used across all metric modules for nucleotide validation

### `metrics/basic.rs` - Basic Metrics
- `snps_distance()` → count (u32 as f64)
- `hamming_distance()` → proportion [0, 1]
- `p_distance()` → proportion [0, 1]
- All skip gaps (`-`) and `N` positions
- All use byte-level iteration for performance

### `metrics/evolutionary.rs` - Evolutionary Models
- `tn93_distance()` → Tamura-Nei 93 (most sophisticated)
- `jc69_distance()` → Jukes-Cantor
- `k2p_distance()` → Kimura 2-parameter
- `tn92_distance()` → Earlier TN variant
- All use frequency counting + mathematical transformations
- Return distance values accounting for multiple substitutions

### `metrics/diagnostics.rs` - Diagnostic Metrics
- `transitions_count()` → count of purine/pyrimidine same-class changes
- `transversions_count()` → count of cross-class changes
- `ts_tv_ratio()` → transitions / transversions (with inf handling)
- Small & isolated - good for quality checks

### `pairwise.rs` - Batch Operations
- `compute_pairwise_distances(nodes: Array, metric: &str) → Array`
- **JS Interop**: Uses `js_sys::{Array, Object, Reflect}` to handle JS objects
- Extracts `{id, seq}` from input, computes all pairs with chosen metric
- Returns array of `{source, target, distance}` objects
- O(n²) complexity for n sequences

## Key Patterns & Conventions

### WASM Bindings
- `#[wasm_bindgen]` attribute on public functions → exposed to JavaScript
- All exported functions must have no generics, lifetimes except `&str`
- Return types: `f64`, `String`, `Array` (from `js_sys`) - primitive JS types
- Parameters: `&str` for sequences (automatically converted from JS strings)

### Gap/Ambiguity Handling
- **Skip strategy**: Continue loop on gaps or `N` (don't count in denominator)
- Enables multi-species alignment support
- Explicit in each metric function

### Performance Optimizations
- Byte-level iteration (`as_bytes()`) instead of char parsing
- Counted accumulation instead of vector allocation
- Single pass through sequences where possible
- `#[profile.release]` uses `opt-level = 3, lto = true`

### Return Types
- Counts: `f64` (for JS interop consistency)
- Proportions/distances: `f64` [0, 1]
- Edge cases return `1.0` (max distance) or `0.0` (identical)

## Build System

### Build Pipeline
```
src/*.rs → wasm-pack → pkg/
  ├── bx_distance_models.js (WASM loader + bindings)
  ├── bx_distance_models_bg.wasm (compiled binary)
  ├── bx_distance_models.d.ts (TypeScript defs)
  └── package.json (for npm publish)
```

### Commands
- `npm run build` → dev build (fast compile, unoptimized WASM)
- `npm run build:release` → optimized build (slower compile, faster runtime)
- `cargo test` → run Rust unit tests

## Dependencies

### Rust
- **wasm-bindgen**: JS ↔ Rust FFI
- **js-sys**: Low-level JS bindings (Array, Object, Reflect)
- **console_error_panic_hook**: Better panic messages in browser (optional feature)

### No Standard Lib Dependencies
- All math done with primitives (no `libm` for basic ops)
- Distance formulas use manual log/exp where needed

## Common Tasks for Agents

### Adding a New Distance Metric
1. Create function in `metrics/new_metric.rs` with `#[wasm_bindgen]` attribute
2. Add to `metrics/mod.rs` with `pub mod` + `pub use`
3. Re-export in `lib.rs` via `pub use metrics::new_metric::*`
4. Add tests in `tests.rs`
5. Rebuild: `npm run build:release`

### Modifying Core Logic
1. Edit function in `src/metrics/*.rs` or `src/utils.rs`
2. Run `cargo test` to validate
3. Rebuild WASM: `npm run build`
4. Test in JavaScript if cross-language behavior changed

### Optimizing Performance
- Profile with `npm run build:release` then benchmark
- Consider vectorization only if WASM SIMD becomes critical
- Byte iteration is already highly optimized
- Avoid allocations in hot paths

### Understanding Test Flow
- `tests.rs` uses `#[cfg(test)]` blocks with `#[test]` attribute
- Can compare Rust results directly
- WASM bindings tested indirectly via JS test suite (if present)

## File I/O & State
- **No file I/O**: Pure computation library
- **No global state**: All functions are stateless, pure
- **Thread-safe**: WASM is single-threaded, but functions have no side effects

## Type Safety Notes
- Rust type system prevents out-of-bounds access
- All string inputs validated via byte iteration
- F64 NaN/Inf handled explicitly in some functions
- JavaScript side should validate `{id, seq}` format before calling `compute_pairwise_distances()`

## Known Limitations
- Sequences compared up to length of shorter sequence
- No support for protein sequences
- No support for RNA (uses DNA alphabet only)
- TN93 distance can be unreliable for very high divergence (p > 0.75)
- Batch operations O(n²) - use streaming for very large sequence sets
