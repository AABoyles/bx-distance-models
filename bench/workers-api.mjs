// Web-Worker-style parallelism with the library's own API: each worker builds a SequenceSet from
// the same strings and computes its share of rows given by row_splits().
//
//   node bench/workers-api.mjs [n=4000] [L=1000] [workerCounts=1 2 4]
//
// Requires `npm run bench:build:api`.
import { Worker, isMainThread, parentPort, workerData } from 'node:worker_threads';
import { createRequire } from 'node:module';
import { makeData } from './data.mjs';

const require = createRequire(import.meta.url);
const bx = require('../pkg-node/bx_distance_models.js');

if (!isMainThread) {
  const { seqs, start, end } = workerData;
  const t0 = performance.now();
  const set = new bx.SequenceSet(seqs);
  const out = set.distances('tn93', start, end);
  set.free();
  parentPort.postMessage({ ms: performance.now() - t0, out }, [out.buffer]);
} else {
  const n = +(process.argv[2] ?? 4000), L = +(process.argv[3] ?? 1000);
  const counts = process.argv.length > 4 ? process.argv.slice(4).map(Number) : [1, 2, 4];
  const seqs = makeData({ n, L });
  const np = n * (n - 1) / 2;

  console.log(`workers (SequenceSet)  n=${n} L=${L} pairs=${np.toLocaleString()}  TN93, full matrix output`);
  for (const k of counts) {
    const set = new bx.SequenceSet(seqs);
    const splits = set.row_splits(k);
    set.free();
    const t0 = performance.now();
    const parts = await Promise.all(Array.from({ length: k }, (_, w) => new Promise((resolve, reject) => {
      const wk = new Worker(new URL(import.meta.url), { workerData: { seqs, start: splits[w], end: splits[w + 1] } });
      wk.once('message', m => { resolve(m); wk.terminate(); });
      wk.once('error', reject);
    })));
    const wall = performance.now() - t0;
    const total = parts.reduce((s, p) => s + p.out.length, 0);
    if (total !== np) throw new Error(`expected ${np} distances, got ${total}`);
    const slowest = Math.max(...parts.map(p => p.ms));
    console.log(`  ${String(k).padStart(2)} worker(s): slowest worker (encode + compute) ${slowest.toFixed(0).padStart(5)} ms, wall incl. spawn + copy-back ${wall.toFixed(0).padStart(5)} ms`);
  }
}
