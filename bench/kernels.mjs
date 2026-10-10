// Kernel benchmark: the library's current byte-loop kernels vs. pre-encoded alternatives.
//
//   node bench/kernels.mjs [--simd] [n=1500] [L=1000] [threshold=0.015]
//
// Requires `npm run bench:build:kernels`. Every alternative is checked against the
// current implementation's output (maxdiff) before its timing is reported.
import { readFileSync } from 'node:fs';
import { makeData } from './data.mjs';

const args = process.argv.slice(2);
const simd = args.includes('--simd');
const [nArg = '1500', LArg = '1000', thrArg = '0.015'] = args.filter(a => a !== '--simd');
const n = +nArg, L = +LArg, thr = +thrArg;
const wasmPath = new URL(simd
  ? './kernels/target/simd/wasm32-unknown-unknown/release/bench_kernels.wasm'
  : './kernels/target/wasm32-unknown-unknown/release/bench_kernels.wasm', import.meta.url);

const { instance } = await WebAssembly.instantiate(readFileSync(wasmPath));
const X = instance.exports;
const mem = () => X.memory.buffer;

const seqs = makeData({ n, L });
const np = n * (n - 1) / 2;
const words = Math.ceil(L / 64);

// Load ASCII once
const pAscii = X.alloc(n * L);
{
  const u8 = new Uint8Array(mem(), pAscii, n * L);
  const enc = new TextEncoder();
  seqs.forEach((s, i) => enc.encodeInto(s, u8.subarray(i * L, (i + 1) * L)));
}
const pOut = X.alloc(np * 8), pOut2 = X.alloc(np * 8), pAll = X.alloc(np * 8 * 3);
const pCodes = X.alloc(n * L), pBits = X.alloc(n * words * 24), pTot = X.alloc(n * 16);

// One warm-up call first: V8 compiles wasm with its baseline tier (Liftoff) and only
// switches a function to the optimizing tier (TurboFan) on a later call. Fast kernels are
// then repeated (~3 s budget, at most 15 runs) and the median is reported, since single
// runs of sub-second kernels vary by +-25% on shared machines.
function time(fn) {
  let t = performance.now();
  fn();
  const warmup = performance.now() - t;
  const reps = Math.min(15, Math.max(1, Math.ceil(3000 / warmup)));
  const runs = [];
  for (let r = 0; r < reps; r++) {
    t = performance.now();
    fn();
    runs.push(performance.now() - t);
  }
  runs.sort((a, b) => a - b);
  return runs[reps >> 1];
}
const f64 = (p, len = np) => new Float64Array(mem(), p, len);
function maxDiff(a, b) {
  let m = 0;
  for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i] - b[i]));
  return m;
}

const rows = [];
const report = (name, ms, baseMs, extra = '') => rows.push({ name, ms, baseMs, extra });

// --- current implementation (verbatim logic from src/metrics)
const tCur = time(() => X.pairs_tn93_current(pAscii, n, L, pOut));
const ref = Float64Array.from(f64(pOut));
report('TN93    current byte loop', tCur, tCur);

const tCurH = time(() => X.pairs_hamming_current(pAscii, n, L, pOut));
const refH = Float64Array.from(f64(pOut));
report('Hamming current byte loop', tCurH, tCurH);

const tCurAll = time(() => X.pairs_all_current(pAscii, n, L, pAll));
const refAll = Float64Array.from(f64(pAll, 3 * np));
report('flat_all current (3 passes/pair)', tCurAll, tCurAll);

// --- pre-encoded bytes + 5x5 count matrix
const tEncL = time(() => X.encode_lut(pAscii, n * L, pCodes));
const tLut = time(() => X.pairs_tn93_lut(pCodes, n, L, pOut2));
report('TN93    LUT 5x5 count matrix', tLut, tCur, `encode ${tEncL.toFixed(1)}ms, maxdiff ${maxDiff(ref, f64(pOut2))}`);

// --- bit-planes
const tEncB = time(() => X.encode_bits_all(pAscii, n, L, pBits, pTot));
const tBits = time(() => X.pairs_tn93_bits(pBits, n, words, pOut2));
report('TN93    bit-planes (10 popcnt/64 sites)', tBits, tCur, `encode ${tEncB.toFixed(1)}ms, maxdiff ${maxDiff(ref, f64(pOut2))}`);
const tFast = time(() => X.pairs_tn93_bits_fast(pBits, pTot, n, words, pOut2));
report('TN93    bit-planes + precomputed freqs', tFast, tCur, `maxdiff ${maxDiff(ref, f64(pOut2))}`);
const tBH = time(() => X.pairs_hamming_bits(pBits, n, words, pOut2));
report('Hamming bit-planes', tBH, tCurH, `maxdiff ${maxDiff(refH, f64(pOut2))}`);
const tBAll = time(() => X.pairs_all_bits(pBits, pTot, n, words, pAll));
report('flat_all bit-planes (1 pass/pair)', tBAll, tCurAll, `maxdiff ${maxDiff(refAll, f64(pAll, 3 * np))}`);

// --- threshold pruning: exact, because TN93 >= p-distance
let kept = 0;
const tThr = time(() => { kept = X.pairs_tn93_bits_rows(pBits, pTot, n, words, 0, n, thr, pOut2); });
const got = f64(pOut2);
let refKept = 0, bad = 0;
for (let i = 0; i < np; i++) {
  if (ref[i] <= thr) { refKept++; if (got[i] !== ref[i]) bad++; }
  else if (got[i] <= thr) bad++;
}
report(`TN93    bit-planes + early exit @${thr}`, tThr, tCur, `kept ${kept}/${np} (ref ${refKept}), mismatches ${bad}`);

console.log(`\n${simd ? 'simd128' : 'baseline'} wasm  n=${n} L=${L} pairs=${np.toLocaleString()}`);
console.log(`${'kernel'.padEnd(42)} ${'ms'.padStart(9)}  ${'ns/pair'.padStart(9)}  ${'ns/site'.padStart(8)}  speedup`);
for (const r of rows) {
  const nsPair = (r.ms * 1e6) / np;
  console.log(`${r.name.padEnd(42)} ${r.ms.toFixed(1).padStart(9)}  ${nsPair.toFixed(1).padStart(9)}  ${(nsPair / L).toFixed(3).padStart(8)}  ${('x' + (r.baseMs / r.ms).toFixed(1)).padStart(7)}  ${r.extra}`);
}
