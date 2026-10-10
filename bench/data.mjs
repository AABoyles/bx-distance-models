// Deterministic synthetic alignment: subtypes -> clusters -> individuals, ts:tv ~2:1,
// ragged gap ends and scattered Ns (ACGTN- alphabet only).
function mulberry32(a) {
  return () => {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const TS = { A: 'G', G: 'A', C: 'T', T: 'C' };
const TV = { A: 'CT', G: 'CT', C: 'AG', T: 'AG' };

export function makeData({ n, L, seed = 42, subtypes = 8, clusters = 12,
  subRate = 0.06, clusterRate = 0.02, indivRate = 0.004, nRate = 0.003, maxRagged = 0.02 }) {
  const rnd = mulberry32(seed);
  const mutate = (s, rate) => s.map(b => {
    if (rnd() >= rate) return b;
    return rnd() < 2 / 3 ? TS[b] : TV[b][(rnd() * 2) | 0];
  });
  const root = Array.from({ length: L }, () => 'ACGT'[(rnd() * 4) | 0]);
  const subs = Array.from({ length: subtypes }, () => mutate(root, subRate));
  const clus = subs.flatMap(s => Array.from({ length: clusters }, () => mutate(s, clusterRate)));
  const seqs = [];
  for (let i = 0; i < n; i++) {
    const s = mutate(clus[(rnd() * clus.length) | 0], indivRate);
    const lead = (rnd() * maxRagged * L) | 0, trail = (rnd() * maxRagged * L) | 0;
    for (let k = 0; k < L; k++) {
      if (k < lead || k >= L - trail) s[k] = '-';
      else if (rnd() < nRate) s[k] = 'N';
    }
    seqs.push(s.join(''));
  }
  return seqs;
}
