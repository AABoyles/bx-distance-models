# BX Distance Models

High-performance nucleotide sequence distance calculations compiled to WebAssembly. This library provides a collection of evolutionary and statistical distance metrics for comparing DNA sequences, optimized for speed and accuracy in browser and Node.js environments.

## Features

- **Multiple Distance Metrics**: Support for basic, evolutionary, and diagnostic distance models
- **WebAssembly Performance**: Compiled to WASM for near-native speed in JavaScript environments
- **Zero Dependencies**: Minimal runtime dependencies for maximum portability
- **Batch Processing**: Efficient pairwise distance computation for multiple sequences
- **Transition/Transversion Analysis**: Built-in diagnostics for evolutionary analysis
- **Gap Handling**: Intelligent handling of gaps and ambiguous nucleotides (N)

## Supported Metrics

### Basic Distance Metrics
- **SNPs Distance**: Count of single nucleotide polymorphisms (absolute count)
- **Hamming Distance**: Proportion of differing positions between sequences (0-1)
- **P-Distance**: Simple proportion of nucleotide differences

### Evolutionary Models
- **TN93 (Tamura-Nei 93)**: Sophisticated model accounting for different transition/transversion rates
- **JC69 (Jukes-Cantor)**: Simple evolutionary distance model
- **K2P (Kimura 2-Parameter)**: Two-parameter model for transitions and transversions
- **TN92 (Tamura-Nei 92)**: Earlier variant of the TN93 model

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

Both sequences must be properly aligned. The library automatically ignores positions with gaps (`-`) or ambiguous nucleotides (`N`):

```javascript
const seq1 = 'ATCG-ATCG';  // Gap at position 4
const seq2 = 'ATCGTNCGG';  // Ambiguous nucleotide (N) at position 5

// Distances are computed ignoring these positions
const distance = hamming_distance(seq1, seq2);
```

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

For computing distances between multiple sequences efficiently:

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
//   { source: 'seq2', target: 'seq3', distance: 1.0 }
// ]

// Using TN93 evolutionary model
const evolutionaryDistances = compute_pairwise_distances(sequences, 'tn93');
```

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

#### `tn93_distance(seq1: string, seq2: string) -> number`
Tamura-Nei 93 distance model (most comprehensive two-parameter model).

#### `jc69_distance(seq1: string, seq2: string) -> number`
Jukes-Cantor 69 distance (simple model, assumes equal substitution rates).

#### `k2p_distance(seq1: string, seq2: string) -> number`
Kimura 2-parameter distance (accounts for different transition/transversion rates).

#### `tn92_distance(seq1: string, seq2: string) -> number`
Earlier Tamura-Nei variant for nucleotide sequences.

### Diagnostics

#### `transitions_count(seq1: string, seq2: string) -> number`
Count of purine-to-purine (A↔G) and pyrimidine-to-pyrimidine (C↔T) substitutions.

#### `transversions_count(seq1: string, seq2: string) -> number`
Count of purine-to-pyrimidine and pyrimidine-to-purine substitutions.

#### `ts_tv_ratio(seq1: string, seq2: string) -> number`
Ratio of transitions to transversions (useful for quality assessment).

### Batch Processing

#### `compute_pairwise_distances(sequences: Array<{id: string, seq: string}>, metric: string) -> Array<{source: string, target: string, distance: number}>`
Efficiently computes all pairwise distances for multiple sequences.

**Parameters:**
- `sequences`: Array of objects with `id` (sequence identifier) and `seq` (DNA sequence)
- `metric`: Distance metric to use ("hamming" or "tn93")

**Returns:** Array of edge objects representing pairwise distances.

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

Run the test suite:
```bash
npm test
```

Or use Rust's native testing:
```bash
cargo test
```

### Project Structure

```
bx-distance-models/
├── src/
│   ├── lib.rs                 # Main library entry point
│   ├── utils.rs               # Utility functions (base validation, etc.)
│   ├── pairwise.rs            # Batch distance computation
│   ├── metrics/
│   │   ├── mod.rs             # Metrics module declarations
│   │   ├── basic.rs           # Basic distance metrics
│   │   ├── evolutionary.rs    # Evolutionary distance models
│   │   └── diagnostics.rs     # Transition/transversion analysis
│   └── tests.rs               # Test suite
├── pkg/                       # Generated WASM bindings
├── Cargo.toml                 # Rust project manifest
├── package.json               # Node.js project manifest
└── README.md                  # This file
```

## Performance Considerations

- **Sequence Length**: Performance scales linearly with sequence length
- **Batch Processing**: Use `compute_pairwise_distances()` for multiple sequences rather than calling single-distance functions repeatedly
- **Memory**: WASM runs in the same heap as JavaScript; large sequences may impact garbage collection
- **Optimized Release Build**: Always use `npm run build:release` for production deployments

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
