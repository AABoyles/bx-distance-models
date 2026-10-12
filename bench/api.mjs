// End-to-end benchmark of the library's public wasm-bindgen API, including JS<->WASM marshalling.
//
//   node bench/api.mjs [n=600] [L=1000] [--legacy-only | --set-only]
//
// Requires `npm run bench:build:api` (a Node-targeted build of the library in pkg-node/).
// "cold" is the first call (may run in V8's baseline tier), "warm" the second.
import { createRequire } from 'node:module';
import { makeData } from './data.mjs';

const require = createRequire(import.meta.url);
const bx = require('../pkg-node/bx_distance_models.js');

const flags = process.argv.slice(2).filter(a => a.startsWith('--'));
const [nArg = '600', LArg = '1000'] = process.argv.slice(2).filter(a => !a.startsWith('--'));
const n = +nArg, L = +LArg;
const seqs = makeData({ n, L });
const nodes = seqs.map((seq, i) => ({ id: `seq${i}`, seq }));
const ids = nodes.map(d => d.id);
const np = n * (n - 1) / 2;

function timed(fn) {
  const t0 = performance.now();
  const result = fn();
  return [performance.now() - t0, result];
}
function run(label, fn) {
  const [cold, result] = timed(fn);
  const [warm] = timed(fn);
  console.log(`${label.padEnd(52)} cold ${cold.toFixed(0).padStart(7)} ms   warm ${warm.toFixed(0).padStart(7)} ms   ${(warm * 1e6 / np).toFixed(0).padStart(6)} ns/pair`);
  return result;
}

console.log(`public API  n=${n} L=${L} pairs=${np.toLocaleString()}`);

if (!flags.includes('--set-only')) {
  for (const metric of ['hamming', 'tn93']) {
    const single = bx[`${metric}_distance`];
    run(`compute_pairwise_distances(objs, '${metric}')`, () => bx.compute_pairwise_distances(nodes, metric));
    run(`compute_pairwise_distances_flat(.., '${metric}')`, () => bx.compute_pairwise_distances_flat(ids, seqs, metric));
    run(`JS loop over ${metric}_distance(a, b)`, () => {
      const out = new Float64Array(np);
      let k = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) out[k++] = single(seqs[i], seqs[j]);
      return out;
    });
  }
  run('compute_pairwise_distances_flat_all(..)', () => bx.compute_pairwise_distances_flat_all(ids, seqs));
}

if (!flags.includes('--legacy-only') && bx.SequenceSet) {
  const [encodeMs, set] = timed(() => new bx.SequenceSet(seqs));
  console.log(`${'new SequenceSet(seqs)  (encode once)'.padEnd(52)} ${encodeMs.toFixed(1)} ms`);
  for (const metric of ['hamming', 'snps', 'tn93']) {
    run(`set.distances('${metric}')`, () => set.distances(metric));
  }
  run(`set.distances('tn93') in blocks of 64 rows`, () => {
    const parts = [];
    for (let r = 0; r < n; r += 64) parts.push(set.distances('tn93', r, Math.min(r + 64, n)));
    return parts;
  });
  for (const [metric, threshold] of [['tn93', 0.015], ['tn93', 0.05], ['hamming', 0.015], ['snps', 15]]) {
    const links = run(`set.links('${metric}', ${threshold})`, () => set.links(metric, threshold));
    console.log(`${''.padEnd(52)} ${links.length.toLocaleString()} links`);
    links.free();
  }
  set.free();

  // Ambiguity fallback: 0.5% of sites are R/Y mixtures in every sequence
  const ambSeqs = makeData({ n, L, ambRate: 0.005 });
  for (const mode of ['skip', 'resolve', 'average']) {
    const ambSet = new bx.SequenceSet(ambSeqs, mode);
    run(`set.distances('tn93'), 0.5% R/Y, ambiguity '${mode}'`, () => ambSet.distances('tn93'));
    ambSet.free();
  }
}
