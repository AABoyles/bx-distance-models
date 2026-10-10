# Benchmarks

Harness behind [`docs/performance-analysis.md`](../docs/performance-analysis.md).
All timings run in Node's V8, the same WebAssembly engine as Chrome.

| Script | What it measures | Build first |
|---|---|---|
| `kernels.mjs` | Per-pair kernels in isolation: verbatim copies of the library's current byte loops vs. pre-encoded alternatives (5×5 count table, 2-bit planes + popcount, early exit). Every alternative is checked against the current output before it is timed. | `npm run bench:build:kernels` |
| `api.mjs` | The library's public wasm-bindgen API end to end, including JS↔WASM marshalling. | `npm run bench:build:api` (needs `wasm-pack`) |
| `workers.mjs` | Splitting the pair triangle across `worker_threads` (same model as Web Workers). | `npm run bench:build:kernels` |
| `data.mjs` | Deterministic synthetic alignment (subtypes → clusters → individuals, ts:tv ≈ 2:1, ragged gap ends, scattered `N`). | — |

```bash
npm run bench:build
npm run bench:kernels               # node bench/kernels.mjs [--simd] [n] [L] [threshold]
npm run bench:api                   # node bench/api.mjs [n] [L]
npm run bench:workers               # node bench/workers.mjs [n] [L] [workerCounts...]
```

`kernels/` is a standalone crate that exports raw `extern "C"` functions (no wasm-bindgen), so the
kernel numbers exclude glue overhead. `bench:build:kernels` builds it twice: baseline and with
`+simd128` (use `--simd` to run the latter). The build scripts use POSIX shell syntax.
