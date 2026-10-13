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
│   ├── encoding.rs               # IUPAC table, 2-bit-plane encoding, sparse ambiguity codes
│   ├── counts.rs                 # Counts struct + kernels (string, bit-plane, ambiguity fallback)
│   ├── sequence_set.rs           # SequenceSet / Links: encode once, dense or thresholded output
│   ├── pairwise.rs               # Legacy batch functions, now thin wrappers over SequenceSet
│   ├── metrics/
│   │   ├── mod.rs                # Module declarations
│   │   ├── basic.rs              # SNPs, Hamming, P-distance
│   │   ├── evolutionary.rs       # TN93, JC69, K2P, T92 (as tn92)
│   │   └── diagnostics.rs        # Transitions/transversions/ratios
│   ├── tests.rs                  # Unit tests (incl. bit-for-bit check against the old kernels)
│   ├── reference_tests.rs        # Checks against ape::dist.dna (table in tests/reference/)
│   └── proptests.rs              # Property-based tests (proptest): kernels vs per-site reference, invariants
├── tests/
│   ├── node/smoke.mjs            # JS binding checks against a Node build (npm run test:node)
│   └── reference/                # generate.R (R + ape) → ape_distances.tsv, read by reference_tests.rs
├── .github/workflows/            # ci.yml (clippy, tests, wasm build, Node checks, bench drivers);
│                                 # reference.yml (regenerates the ape table when it changes)
├── bench/                        # Benchmark harness (see bench/README.md)
│   ├── kernels/                  # Standalone kernel-benchmark crate (raw extern "C" exports)
│   └── *.mjs                     # Node drivers: kernels, public API, workers; synthetic data
├── docs/
│   └── performance-analysis.md   # Where the time goes and what to optimize, with measurements
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

### `encoding.rs` - Sequence Encodings
- `IUPAC[byte]` → 4-bit base set (A=1 C=2 G=4 T=8, N=15, 0 = not a nucleotide); case-insensitive, U = T
- `BASE[byte]` → 2-bit base index (A=0 C=1 G=2 T=3) or 4 for anything else
- `Encoded`: per sequence, `[lo, hi, valid]` bit-planes per 64 sites (A=00 C=01 G=10 T=11, low bit = pyrimidine)
- `AmbiguityCodes`: only for sequences with partial ambiguity codes, and only when the set tracks them:
  a mask plus a sorted list of (position, IUPAC code)
- `Ambiguity`: `Skip` (default), `Resolve`, `Average`

### `counts.rs` - Counting Kernels
- `Counts { sites, ts_purine, ts_pyrimidine, transversions, freq: [A, C, G, T] }` (f64; fractional with ambiguity)
- `Mismatches { sites, mismatches }`: all that SNPs, Hamming, p-distance and JC69 need
- `CountsWithGc { counts, gc1 }`: `Counts` plus the first sequence's G+C count, which only T92 needs.
  Keep `Counts` itself at 64 bytes (8 × f64): a larger struct is copied with `memory.copy`, a
  runtime call in V8, and an extra field cost every full-count metric 5-9%
- **Every metric is a function of `Counts`** (the light ones also of `Mismatches`), so each pair is
  scanned once whatever the metric
- `Tally` trait (`Mismatches`, `Counts`, `CountsWithGc`): `count_planes()` bit-plane kernel, 2 /
  ~10 / ~11 `popcount`s per 64 sites; optional early exit on a mismatch bound
- `Counts::from_strs()` / `CountsWithGc::from_strs()` / `Mismatches::from_strs()` → string kernels
  for the single-pair functions (no branches on bases)
- `Resolutions` + `add_ambiguous_sites()` → ambiguity fallback, visiting only the sparse ambiguous sites
- `count_pair::<T>()` → fast kernel, plus the fallback when either sequence carries ambiguity codes

### `metrics/*.rs` - Metric Formulas
- Each metric has a `pub(crate) fn name(c: &Counts) -> f64` formula (`&Mismatches` for SNPs,
  Hamming/p-distance and JC69, `&CountsWithGc` for T92) and a `#[wasm_bindgen]` string wrapper
  (`name_distance(seq1, seq2)` etc.) that calls it via the matching `from_strs`
- `basic.rs`: SNPs (count), Hamming / p-distance (identical, proportion [0, 1])
- `evolutionary.rs`: TN93, JC69, K2P, T92 (Tamura 1992, exposed as `tn92`). Uncapped
  substitutions per site; `f64::INFINITY` when a log argument would be non-positive, or is within
  rounding of zero (`<= SATURATION`, 1e-12: exact boundary pairs otherwise gave ~19). Zero-weight
  terms of degenerate base compositions drop out (`ratio`/`weighted_ln`) instead of becoming NaN
- `diagnostics.rs`: transitions, transversions, Ts/Tv ratio
- `mod.rs`: `Metric` enum → `parse()` (names accepted from JS), `scan()` (the cheapest `Scan` that
  serves it: `Mismatches` < `Counts` < `CountsWithGc`), `eval()` / `eval_mismatches()` /
  `eval_with_gc()`, and `mismatch_bound()` for early exit (only metrics provably >= p-distance,
  plus SNPs). `eval()` panics for T92, which must go through `eval_with_gc()`

### `sequence_set.rs` - Batch API
- `SequenceSet::new(seqs, ambiguity?)` encodes once; `distances()` (dense, typed array), `links()`
  (thresholded, with early exit), `distance()`, `row_splits()`
- All batch methods take an optional row range; pair order is the row-major upper triangle
- Dense output is written into the JS `Float64Array` row by row (no second full-size copy)
- Native-testable internals (`from_strs`, `counts`, `pair`, `dense`, `sparse`) are separate from
  the `#[wasm_bindgen]` wrappers, which need a JS host

### `pairwise.rs` - Legacy Batch Functions
- `compute_pairwise_distances()` (objects), `compute_pairwise_distances_flat()`,
  `compute_pairwise_distances_flat_all()`, `get_distance_index()`
- Same signatures as before; now built on `SequenceSet` (skip mode). Unknown metric names fall
  back to Hamming

## Key Patterns & Conventions

### WASM Bindings
- `#[wasm_bindgen]` attribute on public functions → exposed to JavaScript
- All exported functions must have no generics, lifetimes except `&str`
- Return types: `f64`, `String`, `Array` (from `js_sys`) - primitive JS types
- Parameters: `&str` for sequences (automatically converted from JS strings)

### Gap/Ambiguity Handling
- One rule for every metric: compare only sites where both sequences hold A/C/G/T (any case,
  U = T); everything else is skipped and not counted in the denominator
- Partial IUPAC codes are skipped by default; `SequenceSet` with `"resolve"`/`"average"` resolves them
- `N` is always skipped

### Performance Optimizations
- Encode once, compare many times: 3 bits per site, 64 sites per `u64` word
- No data-dependent branches in the hot loops (bases are random, so branches on them mispredict)
- One scan per pair for all metrics, and a lighter scan when every requested metric needs only
  mismatches
- Exact early exit for thresholded output (`Metric::mismatch_bound`)
- `count_planes` is `#[inline(never)]` so V8 can tier it up during one long batch call
- `Metric::eval` / `eval_mismatches` are `#[inline(never)]`: inlined into the row loop, LLVM
  hoisted the formulas' `ln` calls (a speculatable intrinsic) above the metric `match`, so every
  pair paid for every metric's logarithms. Check the dense path's speed after touching them
- No JS objects or string copies per pair in the batch API
- Measured in `docs/performance-analysis.md`; `+simd128` was slower, so it is not enabled
- `#[profile.release]` uses `opt-level = 3, lto = true`

### Return Types
- Counts: `f64` (for JS interop consistency)
- Proportions (Hamming, p-distance): `f64` [0, 1]
- Evolutionary distances (TN93, JC69, K2P, T92): substitutions per site, `f64` >= 0 and not capped;
  `Infinity` when the model saturates
- No compared sites at all: `NaN` for every distance (undefined; never within a `links()` threshold).
  SNPs and the transition/transversion counts give `0.0`

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
- `cargo test` (or `npm test`) → unit, ape-reference and property-based tests
- `PROPTEST_CASES=10000 cargo test --release proptests` → longer property-test search
- `cargo clippy --all-targets -- -D warnings` → lint as CI does (also `--target wasm32-unknown-unknown`)
- `Rscript tests/reference/generate.R` → regenerate the ape reference table (needs R + `ape`)
- `npm run test:node` → build for Node and check the JS bindings
- `npm run bench:build` then `npm run bench:kernels` / `bench:api` / `bench:workers` → benchmarks (see `bench/README.md`)

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
1. Write `pub(crate) fn new_metric(c: &Counts) -> f64` in the fitting `metrics/*.rs` file (or take
   `&Mismatches` if compared sites and mismatches are all it needs), plus a `#[wasm_bindgen]`
   string wrapper calling it with `Counts::from_strs` / `Mismatches::from_strs`
2. Add a `Metric` variant in `metrics/mod.rs`: `parse()` name(s), `scan()`, `eval()` (or
   `eval_mismatches()` / `eval_with_gc()` for the other scan levels), `METRIC_NAMES`; and a
   `mismatch_bound()` arm **only** if the metric is provably never below the p-distance (otherwise
   pruning would drop valid links)
3. If it needs more than `Counts` provides, add a new `Tally` type (like `CountsWithGc`) and `Scan`
   level rather than growing `Counts` (see the 64-byte note above), and measure the dense path
4. Add tests in `tests.rs`; `ALL` in the tests (and in `proptests.rs`) lists every metric. If
   `ape::dist.dna` has the model, add its column to `MODELS` in `tests/reference/generate.R`,
   regenerate with `Rscript tests/reference/generate.R`, and map it in `reference_tests.rs`
   (`F81` and `logdet` are already in the table, unused)
5. Rebuild: `npm run build:release`

### Modifying Core Logic
1. Edit the formula in `src/metrics/*.rs` or the kernels in `src/counts.rs`
2. Run `cargo test` to validate
3. Rebuild WASM: `npm run build`
4. Test in JavaScript if cross-language behavior changed

### Optimizing Performance
- Measure with `npm run bench:build` then `npm run bench:api` (and `bench:kernels` for kernel experiments)
- The bit-plane kernel is no longer the bottleneck for typical use; output size and pair count are
- Avoid allocations, JS calls and data-dependent branches in hot paths

### Understanding Test Flow
- `tests.rs`: hand-computed cases, the legacy (pre-encoding) implementation as a reference, and
  fixed-seed random sets
- `reference_tests.rs`: every pair in `tests/reference/ape_distances.tsv`, through both the string
  functions and `SequenceSet`. NA cells are cases ape scores differently by design (see the R
  script's header); saturated models are NaN in ape and `Infinity` here
- `proptests.rs`: generated sets of related sequences in all three ambiguity modes; failures shrink
  to a minimal input. With integer counts (skip mode) checks are bit for bit; with fractional
  counts (resolve/average) they allow rounding, since summation order differs between scan levels
- `tests/node/smoke.mjs`: the wasm-bindgen surface (`npm run test:node`)
- CI (`.github/workflows/ci.yml`) runs all of the above plus clippy and the benchmark drivers

## File I/O & State
- **No file I/O**: Pure computation library
- **No global state**: All functions are stateless, pure; `SequenceSet` is immutable after construction
- **Thread-safe**: WASM is single-threaded, but functions have no side effects

## Type Safety Notes
- Rust type system prevents out-of-bounds access
- All string inputs validated via byte iteration
- F64 NaN/Inf handled explicitly in some functions
- JavaScript side should validate `{id, seq}` format before calling `compute_pairwise_distances()`
- wasm-bindgen methods returning `Result<_, JsError>` throw in JS; keep the fallible logic in
  plain Rust methods (returning `String` errors) so it stays testable with `cargo test`

## Known Limitations
- Sequences compared up to length of shorter sequence
- No support for protein sequences
- RNA is read as DNA (U = T); no protein support
- Evolutionary distances grow quickly near saturation and become `Infinity` beyond it
- Batch operations O(n²) - use `links()` with a threshold, row ranges, or Web Workers for large sets
- A dense matrix must fit in 32-bit WASM memory (about 500 M distances); `distances()` errors beyond that
