# BX Distance Models

High-performance nucleotide sequence distance calculations compiled to WebAssembly. This library provides a collection of evolutionary and statistical distance metrics for comparing DNA sequences, optimized for speed and accuracy in browser and Node.js environments.

## Features

- **Multiple Distance Metrics**: Support for basic, evolutionary, and diagnostic distance models
- **WebAssembly Performance**: Compiled to WASM for near-native speed in JavaScript environments
- **Zero Dependencies**: Minimal runtime dependencies for maximum portability
- **Batch Processing**: `SequenceSet` encodes sequences once (2 bits per base plus a validity mask) and compares 64 sites per instruction with `popcount`
- **Thresholded Output**: return only pairs within a distance threshold, with exact early exit for far-apart pairs
- **Ambiguity Codes**: skip IUPAC ambiguity codes (fast path), or resolve/average them (fallback used only for sequences that contain them)
- **Transition/Transversion Analysis**: Built-in diagnostics for evolutionary analysis
- **Gap Handling**: Gaps, `N` and other non-nucleotide characters are skipped

## Supported Metrics

### Basic Distance Metrics
- **SNPs Distance**: Count of single nucleotide polymorphisms (absolute count)
- **Hamming Distance**: Proportion of differing positions between sequences (0-1)
- **P-Distance**: Simple proportion of nucleotide differences

### Evolutionary Models
- **TN93 (Tamura-Nei 93)**: Sophisticated model accounting for different transition/transversion rates
- **JC69 (Jukes-Cantor)**: Simple evolutionary distance model
- **K2P (Kimura 2-Parameter)**: Two-parameter model for transitions and transversions
- **T92 (Tamura 1992)**, as `tn92`: Kimura's two-parameter model extended for G+C-content bias

Evolutionary distances are expected substitutions per site, so they can exceed 1. When two
sequences are too divergent for a model (one of its logarithms would need a non-positive
argument, e.g. p ≥ 3/4 for JC69), the distance is `Infinity`; it never falls within a finite
`links()` threshold.

### Diagnostics
- **Transitions Count**: Purine-to-purine (A↔G) and pyrimidine-to-pyrimidine (C↔T) changes
- **Transversions Count**: Purine-to-pyrimidine and pyrimidine-to-purine changes
- **Ts/Tv Ratio**: Transition-to-transversion ratio

## Installation

### From NPM
```bash
npm install bx-distance-models
```

### From Source
Clone the repository and follow the [Build for Development](#build-for-development) section.

## Usage

### Basic Example (JavaScript/TypeScript)

```javascript
import { snps_distance, hamming_distance, tn93_distance } from 'bx-distance-models';

// Simple pairwise distances
const seq1 = 'ATCGATCG';
const seq2 = 'ATCGATCC';

// Count SNPs (absolute differences)
const snps = snps_distance(seq1, seq2); // Returns: 1

// Hamming distance (proportion)
const hamming = hamming_distance(seq1, seq2); // Returns: 0.125 (1/8)

// TN93 evolutionary distance
const tn93 = tn93_distance(seq1, seq2); // Returns: evolutionary distance
```

### Handling Gaps and Ambiguous Bases

Both sequences must be properly aligned. Every metric compares only sites where both sequences hold
`A`, `C`, `G` or `T` (either case; `U` is read as `T`). Gaps (`-`, `.`), `N`, `?` and any other
character are skipped and do not count toward the number of compared sites:

```javascript
const seq1 = 'ATCG-ATCG';  // Gap at position 4
const seq2 = 'ATCGTNCGG';  // Ambiguous nucleotide (N) at position 5

// Distances are computed ignoring these positions
const distance = hamming_distance(seq1, seq2);
```

Partial IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) are skipped by
the single-pair functions. `SequenceSet` can resolve them instead; see
[Ambiguity codes](#ambiguity-codes).

### Transition/Transversion Analysis

```javascript
import { transitions_count, transversions_count, ts_tv_ratio } from 'bx-distance-models';

const seq1 = 'ATCG';
const seq2 = 'AGCG';  // A->G transition at position 1

const transitions = transitions_count(seq1, seq2);    // Returns: 1
const transversions = transversions_count(seq1, seq2); // Returns: 0
const ratio = ts_tv_ratio(seq1, seq2);                // Returns: Infinity (1/0)
```

### Batch Processing - Pairwise Distances

For many sequences, build a `SequenceSet`. It encodes every sequence once and returns typed arrays
instead of one JS object per pair:

```javascript
import { SequenceSet, get_distance_index } from 'bx-distance-models';

const seqs = ['ATCGATCG', 'ATCGATCC', 'GCGCGCGC'];
const set = new SequenceSet(seqs);

// Upper triangle, row by row: (0,1), (0,2), (1,2)
const d = set.distances('tn93');            // Float64Array [d01, d02, d12]
get_distance_index(0, 2, seqs.length);      // 1: position of pair (0, 2) in d

// Only pairs within a threshold, as parallel typed arrays.
// For snps/hamming/p_distance/jc69/k2p/tn93, far-apart pairs stop early.
const links = set.links('tn93', 0.015);
const { source, target, distance } = links; // Uint32Array, Uint32Array, Float64Array
links.free();

set.free(); // release WASM memory when done
```

Metric names: `snps`, `hamming`, `p_distance`, `jc69`, `k2p`, `tn92`, `tn93`, `transitions`,
`transversions`, `ts_tv_ratio`.

#### Row ranges and Web Workers

`distances()` and `links()` take an optional row range `[rowStart, rowEnd)`. Results for
consecutive ranges concatenate to the full triangle, so long jobs can be chunked (for progress
or cancellation) or split across Web Workers. `row_splits(k)` gives `k + 1` row boundaries with
nearly equal pair counts:

```javascript
// main.js
const splits = new SequenceSet(seqs).row_splits(navigator.hardwareConcurrency);
const parts = await Promise.all([...Array(splits.length - 1)].map((_, w) => new Promise(resolve => {
  const worker = new Worker(new URL('./distance-worker.js', import.meta.url), { type: 'module' });
  worker.onmessage = e => { resolve(e.data); worker.terminate(); };
  worker.postMessage({ seqs, start: splits[w], end: splits[w + 1] });
})));
// parts[0], parts[1], ... concatenated = the full upper triangle

// distance-worker.js
import { SequenceSet } from 'bx-distance-models';
self.onmessage = ({ data: { seqs, start, end } }) => {
  const set = new SequenceSet(seqs);
  const d = set.distances('tn93', start, end);
  set.free();
  self.postMessage(d, [d.buffer]);
};
```

Each worker encodes its own copy of the sequences (a few milliseconds), so no `SharedArrayBuffer`
or cross-origin isolation is needed.

#### Ambiguity codes

```javascript
const set = new SequenceSet(seqs, 'resolve'); // or 'average'; default 'skip'
```

- `'skip'` (default): sites where either sequence holds a partial ambiguity code are not compared.
  Every pair uses the fast bit-plane kernel.
- `'resolve'`: a site counts as a match if the two codes share a base (`R` vs `A`); otherwise it
  is averaged.
- `'average'`: each site contributes every pair of bases the two codes allow, weighted equally
  (`R` vs `A` is half a match and half an A↔G transition). Counts can be fractional.

With `'resolve'` or `'average'`, sequences without ambiguity codes still use the fast path; only
pairs involving a sequence that contains them pay for the slower fallback, and only in proportion
to the number of ambiguous sites. `set.ambiguous_count()` reports how many sequences use it.
`N` is always skipped.

#### Legacy batch functions

The original batch functions keep their signatures and now run on the same kernel:

```javascript
import { compute_pairwise_distances } from 'bx-distance-models';

const sequences = [
  { id: 'seq1', seq: 'ATCGATCG' },
  { id: 'seq2', seq: 'ATCGATCC' },
  { id: 'seq3', seq: 'GCGCGCGC' }
];

// Compute all pairwise distances using Hamming metric (default)
const distances = compute_pairwise_distances(sequences, 'hamming');
// Returns: [
//   { source: 'seq1', target: 'seq2', distance: 0.125 },
//   { source: 'seq1', target: 'seq3', distance: 1.0 },
//   { source: 'seq2', target: 'seq3', distance: 0.875 }
// ]

// Using TN93 evolutionary model
const evolutionaryDistances = compute_pairwise_distances(sequences, 'tn93');
```

Building one JS object per pair costs several times more than computing the distance, so prefer
`SequenceSet` when there are many sequences. `compute_pairwise_distances_flat(ids, seqs, metric)`
and `compute_pairwise_distances_flat_all(ids, seqs)` return `Float64Array`s in the same order as
`SequenceSet.distances()`.

### React Component Example

```javascript
import React, { useState, useEffect } from 'react';
import { hamming_distance, tn93_distance } from 'bx-distance-models';

function SequenceDistanceCalculator() {
  const [seq1, setSeq1] = useState('ATCGATCG');
  const [seq2, setSeq2] = useState('ATCGATCC');
  const [results, setResults] = useState({});

  useEffect(() => {
    const hamming = hamming_distance(seq1, seq2);
    const tn93 = tn93_distance(seq1, seq2);
    setResults({ hamming, tn93 });
  }, [seq1, seq2]);

  return (
    <div>
      <input 
        value={seq1} 
        onChange={(e) => setSeq1(e.target.value)} 
        placeholder="Sequence 1"
      />
      <input 
        value={seq2} 
        onChange={(e) => setSeq2(e.target.value)} 
        placeholder="Sequence 2"
      />
      <div>
        <p>Hamming Distance: {results.hamming?.toFixed(4)}</p>
        <p>TN93 Distance: {results.tn93?.toFixed(4)}</p>
      </div>
    </div>
  );
}
```

## API Reference

### Basic Metrics

#### `snps_distance(seq1: string, seq2: string) -> number`
Returns the absolute count of single nucleotide polymorphisms.

#### `hamming_distance(seq1: string, seq2: string) -> number`
Returns the proportion of differing positions (0.0 to 1.0).

#### `p_distance(seq1: string, seq2: string) -> number`
Returns the simple proportion of nucleotide differences.

### Evolutionary Models

These return substitutions per site (not capped at 1), or `Infinity` when the model saturates.

#### `tn93_distance(seq1: string, seq2: string) -> number`
Tamura-Nei 93 distance (separate rates for A↔G, C↔T and transversions; base frequencies from the
two sequences).

#### `jc69_distance(seq1: string, seq2: string) -> number`
Jukes-Cantor 69 distance (simple model, assumes equal substitution rates). `Infinity` when p ≥ 3/4.

#### `k2p_distance(seq1: string, seq2: string) -> number`
Kimura 2-parameter distance (accounts for different transition/transversion rates).

#### `tn92_distance(seq1: string, seq2: string) -> number`
Tamura (1992) distance, often called T92: d = −h·ln(1 − P/h − Q) − ½(1 − h)·ln(1 − 2Q), with
h = θ₁ + θ₂ − 2θ₁θ₂ from the G+C content θ of each sequence. Equals K2P when both are 50% G+C.

### Diagnostics

#### `transitions_count(seq1: string, seq2: string) -> number`
Count of purine-to-purine (A↔G) and pyrimidine-to-pyrimidine (C↔T) substitutions.

#### `transversions_count(seq1: string, seq2: string) -> number`
Count of purine-to-pyrimidine and pyrimidine-to-purine substitutions.

#### `ts_tv_ratio(seq1: string, seq2: string) -> number`
Ratio of transitions to transversions (useful for quality assessment).

### Batch Processing

#### `new SequenceSet(seqs: string[], ambiguity?: 'skip' | 'resolve' | 'average')`
Encodes aligned sequences once. Throws on an unknown ambiguity mode.

- `length: number`: number of sequences.
- `distance(metric, i, j) -> number`: one pair.
- `distances(metric, rowStart?, rowEnd?) -> Float64Array`: upper-triangle distances for rows
  `[rowStart, rowEnd)` (default all), in pair order.
- `links(metric, threshold, rowStart?, rowEnd?) -> Links`: pairs with distance `<= threshold`.
- `row_splits(parts) -> Uint32Array`: `parts + 1` row boundaries with balanced pair counts.
- `ambiguous_count() -> number`: sequences that use the ambiguity fallback.
- `free()`: releases WASM memory.

Methods throw on an unknown metric or an invalid row range.

#### `Links`
`length`, and `source` / `target` (`Uint32Array`, sequence indices, `source < target`) and
`distance` (`Float64Array`). Each getter returns a copy, so read it once. Call `free()` when done.

#### `get_distance_index(i: number, j: number, n: number) -> number`
Position of pair `(i, j)`, `i < j`, in a full upper-triangle array for `n` sequences.

#### `compute_pairwise_distances(sequences: Array<{id: string, seq: string}>, metric: string) -> Array<{source: string, target: string, distance: number}>`
Computes all pairwise distances and returns one object per pair.

**Parameters:**
- `sequences`: Array of objects with `id` (sequence identifier) and `seq` (DNA sequence)
- `metric`: Any metric name listed above (unknown names fall back to "hamming")

**Returns:** Array of edge objects representing pairwise distances.

#### `compute_pairwise_distances_flat(ids: string[], seqs: string[], metric: string) -> Float64Array`
Same distances as `SequenceSet.distances(metric)`.

#### `compute_pairwise_distances_flat_all(ids: string[], seqs: string[]) -> Float64Array`
Hamming, TN93 and SNP distances for every pair, concatenated in that order (one scan per pair).

## Build for Development

### Prerequisites
- [Rust](https://www.rust-lang.org/tools/install) (latest stable version)
- [Node.js](https://nodejs.org/) (v14 or later)
- [wasm-pack](https://rustwasm.org/docs/wasm-pack/installation/)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/bx-distance-models.git
cd bx-distance-models
```

2. Install dependencies:
```bash
npm install
```

3. Install Rust and wasm-pack (if not already installed):
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install wasm-pack
```

### Build Commands

**Development build** (faster compilation, slower runtime):
```bash
npm run build
```

**Release build** (slower compilation, optimized for speed):
```bash
npm run build:release
```

This generates WebAssembly files in the `pkg/` directory:
- `bx_distance_models.js` - JavaScript bindings
- `bx_distance_models_bg.wasm` - WebAssembly binary
- `bx_distance_models.d.ts` - TypeScript type definitions

### Testing

Run the Rust test suite (includes a bit-for-bit comparison with the previous implementation):
```bash
npm test        # same as: cargo test
```

Check the JavaScript bindings (typed arrays, optional arguments, thrown errors) against a Node
build:
```bash
npm run test:node
```

### Project Structure

```
bx-distance-models/
├── src/
│   ├── lib.rs                 # Main library entry point
│   ├── encoding.rs            # IUPAC table, bit-plane encoding, ambiguity codes
│   ├── counts.rs              # Per-pair counts and the counting kernels
│   ├── sequence_set.rs        # SequenceSet: encode once, dense/thresholded batch output
│   ├── pairwise.rs            # Legacy batch functions
│   ├── metrics/
│   │   ├── mod.rs             # Metrics module declarations
│   │   ├── basic.rs           # Basic distance metrics
│   │   ├── evolutionary.rs    # Evolutionary distance models
│   │   └── diagnostics.rs     # Transition/transversion analysis
│   └── tests.rs               # Test suite
├── tests/node/smoke.mjs       # JS binding checks (npm run test:node)
├── bench/                     # Benchmark harness
├── docs/                      # Performance analysis
├── pkg/                       # Generated WASM bindings
├── Cargo.toml                 # Rust project manifest
├── package.json               # Node.js project manifest
└── README.md                  # This file
```

## Performance Considerations

- **Use `SequenceSet` for many sequences**: it encodes once and avoids per-pair objects and string copies. Calling the single-pair functions in a JS loop copies both sequences into WASM on every call.
- **Use `links()` with a threshold when you only need close pairs**: far-apart pairs stop early, and the output stays small. A full matrix for 10,000 sequences is 50 M distances (400 MB).
- **Split big jobs**: row ranges let you report progress or spread the work across Web Workers (`row_splits`).
- **SNPs, Hamming, p-distance and JC69 are cheapest**: they only need mismatch counts, so `SequenceSet` computes them 2–3× faster than TN93, T92, K2P and the Ts/Tv metrics, which need full counts
- **Sequence Length**: Performance scales linearly with sequence length
- **Memory**: about 3 bits per site per sequence once encoded (sequences that use the ambiguity fallback add 1 bit per site plus 5 bytes per ambiguous site); call `free()` on `SequenceSet` and `Links` when done
- **Optimized Release Build**: Always use `npm run build:release` for production deployments
- **Benchmarks**: [`docs/performance-analysis.md`](docs/performance-analysis.md) measures where the time goes; the harness lives in [`bench/`](bench/README.md)

## Common Issues

### "Module not found" error
Ensure you've built the project with `npm run build` before importing.

### Type definitions missing
Build with `npm run build` to generate `bx_distance_models.d.ts` in the `pkg/` directory.

### Sequences with different lengths
The library handles sequences of different lengths by comparing only up to the length of the shorter sequence.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Open a pull request

## License

[Add your license information here]

## References

- [TN93 Model](https://academic.oup.com/mbe/article/10/3/512/1048118)
- [Jukes-Cantor Model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969))
- [Kimura 2-Parameter Model](https://academic.oup.com/mbe/article/1/1/1/1011031)
- [WebAssembly](https://webassembly.org/)
- [wasm-bindgen](https://rustwasm.org/docs/wasm-bindgen/)
