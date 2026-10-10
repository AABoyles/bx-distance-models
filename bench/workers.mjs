// Web-Worker-style parallelism: split the pair triangle across worker_threads.
//
//   node bench/workers.mjs [n=4000] [L=1000] [workerCounts=1 2 4]
//
// Each worker owns a private wasm instance and a copy of the bit-packed sequences
// (a few MB at most), so no SharedArrayBuffer / cross-origin isolation is needed.
// Requires `npm run bench:build:kernels`.
import { Worker, isMainThread, parentPort, workerData } from 'node:worker_threads';
import { readFileSync } from 'node:fs';
import { makeData } from './data.mjs';

const wasmBytes = readFileSync(new URL('./kernels/target/wasm32-unknown-unknown/release/bench_kernels.wasm', import.meta.url));
const instance = async () => (await WebAssembly.instantiate(wasmBytes)).instance.exports;

if (!isMainThread) {
  const { bits, tots, n, words, i0, i1, thr } = workerData;
  const X = await instance();
  const pB = X.alloc(bits.byteLength), pT = X.alloc(tots.byteLength);
  new Uint8Array(X.memory.buffer, pB, bits.byteLength).set(new Uint8Array(bits));
  new Uint8Array(X.memory.buffer, pT, tots.byteLength).set(new Uint8Array(tots));
  let pairs = 0;
  for (let i = i0; i < i1; i++) pairs += n - 1 - i;
  const pOut = X.alloc(pairs * 8);
  // Row blocks, no warm-up: V8 has no on-stack replacement for wasm, so repeated calls
  // are what let it tier up (and they give natural progress/cancellation points).
  const t0 = performance.now();
  let kept = 0, off = 0;
  for (let r = i0; r < i1; r += 16) {
    const r1 = Math.min(r + 16, i1);
    kept += X.pairs_tn93_bits_rows(pB, pT, n, words, r, r1, thr, pOut + off * 8);
    for (let i = r; i < r1; i++) off += n - 1 - i;
  }
  const ms = performance.now() - t0;
  const out = new Float64Array(X.memory.buffer, pOut, pairs).slice();
  parentPort.postMessage({ ms, kept, out }, [out.buffer]);
} else {
  const n = +(process.argv[2] ?? 4000), L = +(process.argv[3] ?? 1000);
  const counts = process.argv.length > 4 ? process.argv.slice(4).map(Number) : [1, 2, 4];
  const words = Math.ceil(L / 64), np = n * (n - 1) / 2;
  const X = await instance();
  const seqs = makeData({ n, L });
  const pA = X.alloc(n * L), pB = X.alloc(n * words * 24), pT = X.alloc(n * 16);
  const enc = new TextEncoder();
  const u8 = new Uint8Array(X.memory.buffer, pA, n * L);
  seqs.forEach((s, i) => enc.encodeInto(s, u8.subarray(i * L, (i + 1) * L)));
  X.encode_bits_all(pA, n, L, pB, pT);
  const bits = X.memory.buffer.slice(pB, pB + n * words * 24);
  const tots = X.memory.buffer.slice(pT, pT + n * 16);

  // Balanced split of the triangle by pair count, not by row count.
  function splits(k) {
    const b = [0];
    let acc = 0;
    for (let i = 0; i < n; i++) {
      acc += n - 1 - i;
      if (acc >= (np / k) * b.length && b.length < k) b.push(i + 1);
    }
    b.push(n);
    return b;
  }

  console.log(`workers  n=${n} L=${L} pairs=${np.toLocaleString()}  TN93 bit-planes, full matrix output`);
  for (const k of counts) {
    const b = splits(k);
    const t0 = performance.now();
    const res = await Promise.all(Array.from({ length: k }, (_, w) => new Promise((resolve, reject) => {
      const wk = new Worker(new URL(import.meta.url), {
        workerData: { bits, tots, n, words, i0: b[w], i1: b[w + 1], thr: Infinity },
      });
      wk.once('message', m => { resolve(m); wk.terminate(); });
      wk.once('error', reject);
    })));
    const wall = performance.now() - t0;
    const kern = Math.max(...res.map(r => r.ms));
    console.log(`  ${String(k).padStart(2)} worker(s): slowest kernel ${kern.toFixed(0).padStart(5)} ms, wall incl. spawn+instantiate+copy-back ${wall.toFixed(0).padStart(5)} ms`);
  }
}
