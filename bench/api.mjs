// End-to-end benchmark of the library's public wasm-bindgen API, including JS<->WASM marshalling.
//
//   node bench/api.mjs [n=600] [L=1000]
//
// Requires `npm run bench:build:api` (a Node-targeted build of the library in pkg-node/).
// "cold" is the first call (may run in V8's baseline tier), "warm" the second.
import { createRequire } from 'node:module';
import { makeData } from './data.mjs';

const require = createRequire(import.meta.url);
const bx = require('../pkg-node/bx_distance_models.js');

const n = +(process.argv[2] ?? 600), L = +(process.argv[3] ?? 1000);
const seqs = makeData({ n, L });
const nodes = seqs.map((seq, i) => ({ id: `seq${i}`, seq }));
const ids = nodes.map(d => d.id);
const np = n * (n - 1) / 2;

function timed(fn) {
  const t0 = performance.now();
  fn();
  return performance.now() - t0;
}
function run(label, fn) {
  const cold = timed(fn), warm = timed(fn);
  console.log(`${label.padEnd(48)} cold ${cold.toFixed(0).padStart(7)} ms   warm ${warm.toFixed(0).padStart(7)} ms   ${(warm * 1e6 / np).toFixed(0).padStart(6)} ns/pair`);
}

console.log(`public API  n=${n} L=${L} pairs=${np.toLocaleString()}`);
for (const metric of ['hamming', 'tn93']) {
  const single = metric === 'tn93' ? bx.tn93_distance : bx.hamming_distance;
  run(`compute_pairwise_distances(objs, '${metric}')`, () => bx.compute_pairwise_distances(nodes, metric));
  run(`compute_pairwise_distances_flat(.., '${metric}')`, () => bx.compute_pairwise_distances_flat(ids, seqs, metric));
  run(`JS loop over ${metric}_distance(a, b)`, () => {
    const out = new Float64Array(np);
    let k = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) out[k++] = single(seqs[i], seqs[j]);
    return out;
  });
}
