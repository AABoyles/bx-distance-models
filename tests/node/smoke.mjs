// Exercises the wasm-bindgen surface that `cargo test` cannot reach (typed arrays, optional
// arguments, thrown errors). Run with `npm run test:node`.
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';

const require = createRequire(import.meta.url);
const bx = require('../../pkg-node/bx_distance_models.js');

const seqs = ['ACGTACGTAC', 'ACGTTCGTAC', 'ACRTACGTAC', 'TTGTACGAAC', 'acgtacgtac'];
const n = seqs.length;
const pairs = [];
for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) pairs.push([i, j]);

function check(name, fn) {
  fn();
  console.log(`ok - ${name}`);
}

const set = new bx.SequenceSet(seqs);

check('dense distances match the single-pair functions', () => {
  assert.equal(set.length, n);
  for (const metric of ['hamming', 'snps', 'tn93', 'jc69', 'k2p']) {
    const d = set.distances(metric);
    assert.ok(d instanceof Float64Array);
    assert.equal(d.length, pairs.length);
    const single = bx[`${metric}_distance`];
    pairs.forEach(([i, j], k) => {
      assert.equal(d[k], single(seqs[i], seqs[j]), `${metric} ${i},${j}`);
      assert.equal(set.distance(metric, i, j), d[k]);
      assert.equal(bx.get_distance_index(i, j, n), k);
    });
  }
});

check('row ranges from row_splits concatenate to the full triangle', () => {
  const full = Array.from(set.distances('tn93'));
  for (const parts of [1, 2, 3, 8]) {
    const splits = set.row_splits(parts);
    assert.ok(splits instanceof Uint32Array);
    assert.equal(splits.length, parts + 1);
    const joined = [];
    for (let k = 0; k < parts; k++) joined.push(...set.distances('tn93', splits[k], splits[k + 1]));
    assert.deepEqual(joined, full);
  }
});

check('links are the dense pairs within the threshold', () => {
  const dense = set.distances('hamming');
  const links = set.links('hamming', 0.15);
  const [source, target, distance] = [links.source, links.target, links.distance];
  assert.ok(source instanceof Uint32Array && target instanceof Uint32Array && distance instanceof Float64Array);
  const expected = pairs.filter((_, k) => dense[k] <= 0.15);
  assert.equal(links.length, expected.length);
  expected.forEach(([i, j], k) => {
    assert.deepEqual([source[k], target[k]], [i, j]);
    assert.equal(distance[k], dense[bx.get_distance_index(i, j, n)]);
  });
  links.free();
});

check('ambiguity modes', () => {
  // seqs[2] has R (A|G) where seqs[0] has G: skipped by default, a match when resolved
  assert.equal(set.ambiguous_count(), 0);
  const resolve = new bx.SequenceSet(seqs, 'resolve');
  assert.equal(resolve.ambiguous_count(), 1);
  assert.equal(resolve.distance('snps', 0, 2), 0);
  const average = new bx.SequenceSet(seqs, 'average');
  assert.equal(average.distance('snps', 0, 2), 0.5);
  assert.equal(average.distance('transitions', 0, 2), 0.5);
  resolve.free();
  average.free();

  // GapMM: a gap facing a base is an N (3/4 mismatch, 3 such sites); gap facing gap is skipped
  const gapmm = new bx.SequenceSet(['ACGT-', 'AC--A', 'ACGT-'], 'gapmm');
  assert.equal(gapmm.ambiguous_count(), 3);
  assert.equal(gapmm.distance('snps', 0, 1), 2.25);
  assert.equal(gapmm.distance('hamming', 0, 1), 2.25 / 5);
  assert.equal(gapmm.distance('snps', 0, 2), 0);
  assert.ok(gapmm.distance('tn93', 0, 1) > 0);
  gapmm.free();

  // N and ? are any base outside skip mode: a match when resolved, 3/4 mismatch when averaged
  const withN = ['ACGTA', 'ACGNA', 'AC?TA'];
  const nSkip = new bx.SequenceSet(withN);
  const nResolve = new bx.SequenceSet(withN, 'resolve');
  const nAverage = new bx.SequenceSet(withN, 'average');
  assert.equal(nSkip.ambiguous_count(), 0);
  assert.equal(nResolve.ambiguous_count(), 2);
  assert.equal(nSkip.distance('hamming', 0, 1), 0);
  assert.equal(nResolve.distance('hamming', 0, 1), 0);
  assert.equal(nResolve.distance('hamming', 0, 2), 0);
  assert.equal(nAverage.distance('snps', 0, 1), 0.75);
  assert.equal(nAverage.distance('hamming', 0, 2), 0.75 / 5);
  [nSkip, nResolve, nAverage].forEach(s => s.free());

  // tn93 -a resolve -g F: R facing A is 1 resolved site of 10, so averaged at 0.1, resolved above
  const gSeqs = ['ACGTACGTAC', 'RCGTACGTAC'];
  const g = f => new bx.SequenceSet(gSeqs, 'resolve', { maxAmbiguityFraction: f });
  const [low, high] = [g(0.1), g(0.11)];
  assert.equal(low.distance('snps', 0, 1), 0.5);
  assert.equal(high.distance('snps', 0, 1), 0);
  assert.equal(low.links('snps', 0).length, 0);
  assert.equal(high.links('snps', 0).length, 1);
  [low, high].forEach(s => s.free());
  new bx.SequenceSet(gSeqs, 'resolve', {}).free();

  // tn93K2pFallback: a pair with only A and T uses K2P, as tn93 does
  const atSeqs = ['ATTAATATTA', 'ATTAATTTTA'];
  const [plain, k2pFallback] = [new bx.SequenceSet(atSeqs), new bx.SequenceSet(atSeqs, 'skip', { tn93K2pFallback: true })];
  assert.equal(k2pFallback.distance('tn93', 0, 1), plain.distance('k2p', 0, 1));
  assert.notEqual(plain.distance('tn93', 0, 1), plain.distance('k2p', 0, 1));
  assert.deepEqual(Array.from(k2pFallback.links('tn93', 10).distance), [plain.distance('k2p', 0, 1)]);
  [plain, k2pFallback].forEach(s => s.free());
  new bx.SequenceSet(gSeqs, 'resolve', { maxAmbiguityFraction: undefined }).free();
});

check('errors are thrown as JS errors', () => {
  assert.throws(() => new bx.SequenceSet(seqs, 'bogus'), /unknown ambiguity mode/);
  assert.throws(() => new bx.SequenceSet(seqs, 'resolve', { maxAmbiguity: 0.1 }), /unknown option 'maxAmbiguity'/);
  assert.throws(() => new bx.SequenceSet(seqs, 'average', { maxAmbiguityFraction: 0.1 }), /needs ambiguity mode 'resolve'/);
  assert.throws(() => new bx.SequenceSet(seqs, 'resolve', { maxAmbiguityFraction: 2 }), /between 0 and 1/);
  assert.throws(() => new bx.SequenceSet(seqs, 'resolve', { maxAmbiguityFraction: '0.1' }), /must be a number/);
  assert.throws(() => new bx.SequenceSet(seqs, 'skip', { tn93K2pFallback: 1 }), /must be a boolean/);
  assert.throws(() => set.distances('bogus'), /unknown metric/);
  assert.throws(() => set.distances('tn93', 3, 2), /invalid row range/);
  assert.throws(() => set.distance('tn93', 0, n), /out of range/);
});

check('legacy batch functions', () => {
  const nodes = seqs.map((seq, i) => ({ id: `s${i}`, seq }));
  const tn93 = set.distances('tn93');
  const objs = bx.compute_pairwise_distances(nodes, 'tn93');
  assert.equal(objs.length, pairs.length);
  pairs.forEach(([i, j], k) => assert.deepEqual(objs[k], { source: `s${i}`, target: `s${j}`, distance: tn93[k] }));

  assert.deepEqual(Array.from(bx.compute_pairwise_distances_flat([], seqs, 'tn93')), Array.from(tn93));
  const all = bx.compute_pairwise_distances_flat_all([], seqs);
  const np = pairs.length;
  assert.deepEqual(Array.from(all.subarray(0, np)), Array.from(set.distances('hamming')));
  assert.deepEqual(Array.from(all.subarray(np, 2 * np)), Array.from(tn93));
  assert.deepEqual(Array.from(all.subarray(2 * np)), Array.from(set.distances('snps')));
  assert.equal(bx.compute_pairwise_distances_flat([], [], 'tn93').length, 0);
});

set.free();
