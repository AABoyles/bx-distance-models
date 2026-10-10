# Performance analysis: in-browser pairwise distance computation

**Question.** How can in-browser distance computation be made faster? The working hypothesis was a
compact nucleotide encoding with bitwise operations, in the style of
[NtSeq](https://github.com/keithwhor/NtSeq).

**Short answer.** The hypothesis is right about the kernel, and by more than expected: a
2-bit-plane encoding with `popcount` computes TN93 **~80× faster** than the current byte loop, with
bit-identical results. But on its own it would only make the documented batch API
(`compute_pairwise_distances`) about **5× faster**, because that function spends ~3.3 µs per pair
building JS objects, roughly 20× what the bit-packed kernel needs. The output path and the
encoding have to change together. After that, exact threshold pruning and Web Workers multiply the
gain.

All numbers below come from the harness in [`bench/`](../bench/README.md) and can be reproduced
with `npm run bench:build && npm run bench:kernels`.

## Method

- **Engine:** Node 22 (V8, the same WebAssembly engine as Chrome), one core of a 2.8 GHz Xeon,
  4 vCPUs for the worker test. Absolute numbers will differ on other machines; the ratios are what
  matter.
- **Data:** deterministic synthetic alignment ([`bench/data.mjs`](../bench/data.mjs)) shaped like
  HIV *pol* surveillance data: 8 subtypes (~6% apart) → 12 clusters each (~2%) → individuals
  (~0.4%), transition:transversion ≈ 2:1, ragged gap ends (up to 2% per end), 0.3% `N`.
  Alphabet is `ACGTN-` only, so every metric in the library agrees on which sites to compare.
- **Kernels** ([`bench/kernels/src/lib.rs`](../bench/kernels/src/lib.rs)) are exported as raw
  `extern "C"` functions so the timings exclude wasm-bindgen glue. The "current" kernels are
  verbatim copies of `src/metrics/*.rs`. Every alternative is compared against the current
  output: all of them reproduce it exactly (max |difference| = 0).
- Each kernel is called once to warm up before it is timed, because V8 compiles WebAssembly
  with a baseline tier (Liftoff) first and switches to the optimizing tier (TurboFan) only on a
  later call. Kernels that finish in under a second are then repeated and the median is
  reported: single runs of the fast kernels varied by ±25% on this shared machine.
- The library build used for the API numbers is `wasm-pack build --release` (including
  `wasm-opt`), as for a published package.

## Results

### Kernels (n = 1,500, L = 1,000, 1.12 M pairs, one core)

| Kernel | ns/pair | ns/site | Speedup |
|---|---:|---:|---:|
| TN93, current byte loop | 14,810 | 14.8 | 1× |
| TN93, encode once to 1 byte/site + branch-free 5×5 count table | 1,322 | 1.3 | 11× |
| **TN93, 2 bit-planes + validity mask, `popcount`** | **176** | **0.18** | **84×** |
| TN93, bit-planes + precomputed per-sequence base totals | 220 | 0.22 | 67× |
| TN93, bit-planes + exact early exit at d > 0.015 | 78 | 0.08 | 190× |
| Hamming, current byte loop | 3,386 | 3.4 | 1× |
| Hamming, bit-planes | 41 | 0.04 | 83× |
| `flat_all` (Hamming + TN93 + SNPs), current: 3 passes per pair | 21,927 | 21.9 | 1× |
| `flat_all`, bit-planes: 1 pass per pair | 268 | 0.27 | 82× |

Encoding is negligible: ~5 ms for all 1,500 sequences, against seconds of pairwise work.

**Long sequences** (n = 250, L = 10,000): bit-plane TN93 is 149× faster than the current loop,
and 570× with early exit. The current loop stays at ~15 ns/site, so its cost is per-site work, not
cache misses.

### The public API (n = 600, L = 1,000, 180 K pairs, warm call)

| Call | ns/pair |
|---|---:|
| `compute_pairwise_distances(objs, 'tn93')` | 18,752 |
| `compute_pairwise_distances_flat(ids, seqs, 'tn93')` | 14,705 |
| JS loop over `tn93_distance(a, b)` | 22,897 |
| `compute_pairwise_distances(objs, 'hamming')` | 6,690 |
| `compute_pairwise_distances_flat(ids, seqs, 'hamming')` | 3,495 |
| JS loop over `hamming_distance(a, b)` | 12,139 |

With very short sequences (n = 1,500, L = 16), where the kernel costs almost nothing, the object API
still takes **3,323 ns/pair** (Hamming) against **78 ns/pair** for the flat one. That ~3.3 µs per pair
is pure marshalling: one JS object per pair, three `Reflect.set` calls with freshly converted key
strings, and two id strings converted from Rust to JS for every pair. The README recommends this
function for batch work, and it is the slowest path. The JS loop over the single-pair functions
pays for copying both sequences from JS strings into WASM memory on every call (~8 µs/pair at
L = 1,000).

### Web Workers (n = 4,000, L = 1,000, 8 M pairs, bit-plane TN93)

| Workers | Slowest worker | Wall time incl. spawn, instantiate, copy back |
|---:|---:|---:|
| 1 | 1,998 ms | 2,185 ms |
| 2 | 1,071 ms | 1,184 ms |
| 4 | 497 ms | 610 ms |

3.6× wall-clock on 4 cores (4.0× for the kernels alone). Each worker gets its own WASM instance
and a copy of the encoded sequences (1.5 MB here), so no `SharedArrayBuffer` or cross-origin
isolation headers are needed. The triangle is split by pair count, not by row count. For comparison, the current object API would need about
150 s for these 8 M pairs (extrapolated from 18.75 µs/pair).

## Why the current loop is slow

~15 ns per site is about 40 CPU cycles. The cause is the frequency counting in
`src/metrics/evolutionary.rs`:

```rust
match a { b'A' => freq_a += 1, b'C' => freq_c += 1, b'G' => freq_g += 1, b'T' => freq_t += 1, _ => {} }
match b { /* same */ }
```

Each `match` branches on which base it sees. Bases are effectively random, so the CPU mispredicts
these branches most of the time, twice per site. The evidence: the 1-byte table kernel does the
same byte-at-a-time work but replaces the branches with an indexed increment
(`counts[code1 * 5 + code2] += 1`), and that alone is 11× faster. Hamming has no such `match` and
runs at 3.4 ns/site, still branchy but predictable.

## Encoding design

### NtSeq

NtSeq uses 4 bits per base, one bit per nucleotide (A = `1000`, T = `0100`, G = `0010`,
C = `0001`). An ambiguity code is the OR of its bases (R = A|G, N = all four), and a gap is 0. Eight
bases fit in a 32-bit integer. Matching is a bitwise AND (a nonzero nibble means "compatible") plus
a lookup-table popcount, because JavaScript has no popcount instruction. That design suits
sequence search, where `N` should match anything.

### Recommended for distances: 2 bit-planes + validity mask

```text
A = 00   C = 01   G = 10   T = 11        (low bit = pyrimidine)
per 64 sites:  lo: u64   hi: u64   valid: u64
```

For each 64-site word of a pair:

```text
v            = valid1 & valid2
transversion = (lo1 ^ lo2) & v                      // purine <-> pyrimidine
transition   = (hi1 ^ hi2) & v & ~transversion      // same class, different base
  A<->G      = transition & ~lo1
  C<->T      = transition &  lo1
base freqs   = popcount(v & lo), popcount(v & hi), popcount(v & lo & hi)  for each sequence
```

About 10 `popcount`s per 64 sites, with no data-dependent branches. WebAssembly has a native
`i64.popcnt` instruction (it compiles to x86 `POPCNT`), so the idea works better in WASM than it
did in JavaScript. Storage is 3 bits per site instead of 8: 10,000 sequences of 1 kb take
~3.8 MB, which stays cache-friendly.

**Trade-off: ambiguity codes.** TN93 currently skips any site that is not `A/C/G/T`, so the
2-plane layout reproduces it exactly. HIV-TRACE-style ambiguity resolution (R and Y are common in
HIV population sequencing) needs a richer encoding. NtSeq's 4-bit one-hot code represents it
natively, at more bits per site and more operations per transition test. A good compromise is to
keep the 2-plane fast path for every sequence and store a 4-bit IUPAC code only for sequences
that actually contain ambiguities, resolving those (rare) sites in a separate fix-up pass.

## What did not help

- **`+simd128`.** Compiling the same kernels with WebAssembly SIMD made the bit-plane TN93 kernel
  *slower* (176 → 279 ns/pair). LLVM vectorizes the 64-bit popcount with `i8x16.popcnt` plus
  widening adds, which is worse than the scalar instruction. Hand-written SIMD might win
  something, but the kernel is no longer the bottleneck. Don't start there.
- **Precomputing per-sequence base totals.** This cuts the popcounts per word from 10 to 4, but the
  branch it needs (to correct for sites valid in only one sequence) made the kernel slower
  (220 vs. 176 ns/pair). Popcounts are cheap; branches are not.
- **One giant call.** V8 does not switch a WebAssembly function to optimized code while a call is
  in progress (no on-stack replacement). A single long call therefore runs in baseline code from
  start to finish: the Hamming flat API took 1,117 ms on its first call and 628 ms on its second.
  With baseline code only (`node --liftoff-only`), the bit-plane kernel takes 287 instead of
  176 ns/pair. Long work should be split into calls over blocks of rows, which is also what
  progress reporting and cancellation need.

## Recommendations, in order of payoff

1. **Fix the output path.** Stop building one JS object per pair in Rust. Return typed arrays (a
   `Float32Array` would halve memory), or index pairs plus distances, and build objects in JS only
   if needed. Better still, return only pairs under a threshold: a dense matrix for n = 10,000 is
   50 M pairs = 400 MB as `f64`, and twice that while `compute_pairwise_distances_flat` copies its
   `Vec` into a `Float64Array`. That is a hard limit in 32-bit WebAssembly.
2. **Encode once.** A `SequenceSet` object that takes the sequences once and stores them as
   bit-planes. This also removes the per-call string copies.
3. **One counting kernel for every metric.** Every metric in the library is a cheap function of
   the same eight counts: compared sites, A↔G and C↔T transitions, transversions, and the four
   base frequencies. Counting once replaces `flat_all`'s three passes, `ts_tv_ratio`'s two, and
   the extra pass when TN93 falls back to Hamming.
4. **Exact early exit when a threshold is given.** TN93 is never smaller than the p-distance:
   applying −ln(1−x) ≥ x to each log term gives d ≥ P₁ + P₂ + c·Q, and the coefficient c works
   out to exactly 1, so d ≥ P₁ + P₂ + Q = p. This still holds with the current clamp to 1.0 and
   the fallback to Hamming, with one degenerate exception in the current code (see the last
   section). JC69 and K2P satisfy the same bound. So a pair can stop as soon as its mismatches so
   far exceed `threshold × min(valid sites in either sequence)`. In the benchmark this kept
   exactly the same 11,416 pairs with identical values.
5. **Web Workers.** Split the triangle by pair count; one WASM instance per worker. 3.6× on 4
   cores, wall-clock.
6. **Later, only if profiling points back at the kernel:** hand-written SIMD, cache blocking.

## Behavior decisions an encoding forces

- **Which characters count as data.** Hamming, SNPs, JC69, K2P and the transition/transversion
  counts skip only `-` and `N`; TN93, TN92 and p-distance skip anything that is not upper-case
  `ACGT`. So `hamming_distance("acgt", "ACGT")` is 1.0, `R` vs `A` is a Hamming mismatch but is
  skipped by TN93, and when TN93 falls back to Hamming it mixes the two rules. An encoder has to
  pick one rule (A/C/G/T in either case, U read as T, everything else is not a base), which
  changes Hamming results for inputs containing lower case or ambiguity codes.
- **A TN93 edge case.** If A or G (or C or T) is absent from every compared site, `a1` (or `a2`)
  is 0 and the current code computes 0/0. The NaN then passes through `d.max(0.0)`, which returns
  0.0. Example: `tn93_distance("AACCTTAACCTT", "AACCTTAACCTA")` is 0.0 despite a mismatch. The
  limit of that term is 0, so the fix is to drop it rather than divide.
- **Unrelated to speed, worth checking before locking in tests:** `tn92_distance` does not match the
  Tamura (1992) formula, which is parameterized by GC content. TN93 also caps distances at 1.0 and
  falls back to the p-distance when the model saturates, so the most divergent pairs can report a
  smaller distance than moderately divergent ones.

## Reproducing

```bash
npm run bench:build                     # kernels crate (baseline + simd128) and a Node build of the library
node bench/kernels.mjs 1500 1000        # kernel table
node bench/kernels.mjs --simd 1500 1000 # same kernels compiled with +simd128
node bench/kernels.mjs 250 10000        # long sequences
node bench/api.mjs 600 1000             # public API, realistic length
node bench/api.mjs 1500 16              # public API, marshalling overhead only
node bench/workers.mjs 4000 1000 1 2 4  # worker scaling
node --liftoff-only bench/kernels.mjs 600 1000   # baseline tier only
```
