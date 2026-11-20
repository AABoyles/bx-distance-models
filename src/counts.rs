//! Per-pair site counts and the kernels that produce them.
//!
//! Every metric in this crate is a function of [`Counts`], so each pair is scanned once no
//! matter how many metrics are wanted. Metrics that only need compared sites and mismatches
//! (SNPs, Hamming, p-distance, JC69) can use the cheaper [`Mismatches`] scan instead.
use crate::encoding::{Ambiguity, Encoded, FullWord, Word, ANY, BASE};

/// Base indices: A=0 C=1 G=2 T=3 (purines even, pyrimidines odd).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Counts {
    /// Sites compared (both sequences hold a base).
    pub sites: f64,
    /// A <-> G transitions.
    pub ts_purine: f64,
    /// C <-> T transitions.
    pub ts_pyrimidine: f64,
    pub transversions: f64,
    /// A, C, G, T counts over the compared sites of both sequences (sums to `2 * sites`).
    pub freq: [f64; 4],
}

/// [`Counts`] plus the G+C content of the first sequence, which only T92 needs.
///
/// Kept out of `Counts` on purpose: tracking it in every full scan cost the other full-count
/// metrics 10-25% (one more accumulator in the hot loop), and even an unused field cost 5-9%,
/// because a struct over 64 bytes is copied with `memory.copy`, a runtime call in V8.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CountsWithGc {
    pub counts: Counts,
    /// G or C in the first sequence at compared sites (the second's is `freq[C] + freq[G] - gc1`).
    pub gc1: f64,
}

impl Counts {
    #[inline]
    pub fn transitions(&self) -> f64 {
        self.ts_purine + self.ts_pyrimidine
    }

    #[inline]
    pub fn mismatches(&self) -> f64 {
        self.transitions() + self.transversions
    }

    /// Counts two strings site by site, up to the shorter length, skipping non-bases.
    pub fn from_strs(seq1: &str, seq2: &str) -> Counts {
        CountsWithGc::from_strs(seq1, seq2).counts
    }

    /// Adds `w` sites where sequence 1 holds base `x` and sequence 2 holds base `y`.
    #[inline]
    pub fn add(&mut self, x: usize, y: usize, w: f64) {
        self.sites += w;
        self.freq[x] += w;
        self.freq[y] += w;
        if x != y {
            if (x ^ y) & 1 == 1 {
                self.transversions += w;
            } else if x & 1 == 0 {
                self.ts_purine += w;
            } else {
                self.ts_pyrimidine += w;
            }
        }
    }

    #[inline]
    fn add_counts(&mut self, other: &Counts) {
        self.add_counts_scaled(other, 1.0);
    }

    #[inline]
    fn add_counts_scaled(&mut self, other: &Counts, n: f64) {
        if n == 0.0 {
            return;
        }
        self.sites += n * other.sites;
        self.ts_purine += n * other.ts_purine;
        self.ts_pyrimidine += n * other.ts_pyrimidine;
        self.transversions += n * other.transversions;
        for (f, o) in self.freq.iter_mut().zip(other.freq) {
            *f += n * o;
        }
    }
}

impl CountsWithGc {
    /// Counts two strings site by site, up to the shorter length, skipping non-bases.
    pub fn from_strs(seq1: &str, seq2: &str) -> CountsWithGc {
        // 5x5 matrix indexed by base (A C G T, 4 = not a base): one indexed increment per
        // site, no data-dependent branches.
        let mut m = [0u32; 32];
        for (&a, &b) in seq1.as_bytes().iter().zip(seq2.as_bytes()) {
            m[(BASE[a as usize] as usize * 5 + BASE[b as usize] as usize) & 31] += 1;
        }
        let mut c = CountsWithGc::default();
        for x in 0..4 {
            for y in 0..4 {
                let n = m[x * 5 + y];
                if n != 0 {
                    c.add(x, y, n as f64);
                }
            }
        }
        c
    }

    /// Adds `w` sites where sequence 1 holds base `x` and sequence 2 holds base `y`.
    #[inline]
    pub fn add(&mut self, x: usize, y: usize, w: f64) {
        self.counts.add(x, y, w);
        if matches!(x, 1 | 2) {
            self.gc1 += w;
        }
    }

    /// Adds `n` sites that each contribute `site`.
    #[inline]
    fn add_scaled(&mut self, site: &CountsWithGc, n: f64) {
        self.counts.add_counts_scaled(&site.counts, n);
        self.gc1 += n * site.gc1;
    }

    /// Adds `w` sites holding IUPAC codes `c1` and `c2`, resolving ambiguity per `mode`.
    pub fn add_codes(&mut self, c1: u8, c2: u8, w: f64, mode: Ambiguity) {
        let shared = c1 & c2;
        if mode == Ambiguity::Resolve && shared != 0 {
            let share = w / shared.count_ones() as f64;
            for b in bases(shared) {
                self.add(b, b, share);
            }
        } else {
            let share = w / (c1.count_ones() * c2.count_ones()) as f64;
            for x in bases(c1) {
                for y in bases(c2) {
                    self.add(x, y, share);
                }
            }
        }
    }
}

/// Compared sites and the mismatches among them: all that SNPs, Hamming, p-distance and JC69 need.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Mismatches {
    pub sites: f64,
    pub mismatches: f64,
}

impl Mismatches {
    /// Counts two strings site by site, up to the shorter length, skipping non-bases.
    pub fn from_strs(seq1: &str, seq2: &str) -> Mismatches {
        let (mut sites, mut mismatches) = (0u32, 0u32);
        for (&a, &b) in seq1.as_bytes().iter().zip(seq2.as_bytes()) {
            let (x, y) = (BASE[a as usize], BASE[b as usize]);
            // Bases are 0..=3 and anything else is 4, so both are bases exactly when x | y < 4.
            let compared = ((x | y) < 4) as u32;
            sites += compared;
            mismatches += compared & (x != y) as u32;
        }
        Mismatches { sites: sites as f64, mismatches: mismatches as f64 }
    }
}

impl From<&Counts> for Mismatches {
    fn from(c: &Counts) -> Self {
        Mismatches { sites: c.sites, mismatches: c.mismatches() }
    }
}

fn bases(code: u8) -> impl Iterator<Item = usize> {
    (0..4).filter(move |b| code >> b & 1 == 1)
}

/// What one site contributes for every pair of IUPAC codes under one ambiguity mode, built once
/// per `SequenceSet` so the per-site work in the fallback is a lookup and eight additions.
pub(crate) struct Resolutions(Box<[CountsWithGc; 256]>);

impl Resolutions {
    pub fn new(mode: Ambiguity) -> Self {
        let mut table = Box::new([CountsWithGc::default(); 256]);
        // Code 0 (not a base) is never compared. N (15) is only looked up in modes that compare it.
        for c1 in 1..=ANY {
            for c2 in 1..=ANY {
                table[(c1 as usize) << 4 | c2 as usize].add_codes(c1, c2, 1.0, mode);
            }
        }
        let table = Resolutions(table);
        // `Mismatches::add_full_sites` relies on a base facing an N weighing the same whatever the
        // base, and on either side.
        for b in [1u8, 2, 4, 8] {
            for site in [table.get(b, ANY), table.get(ANY, b)] {
                debug_assert_eq!(Mismatches::from(&site.counts), Mismatches::from(&table.get(1, ANY).counts));
            }
        }
        table
    }

    #[inline]
    fn get(&self, c1: u8, c2: u8) -> &CountsWithGc {
        &self.0[(c1 as usize) << 4 | c2 as usize]
    }
}

/// Mismatches are checked for early exit every 8 words (512 sites).
const CHUNK: usize = 8;

/// What a scan of one pair accumulates: full [`Counts`], or only [`Mismatches`].
pub(crate) trait Tally: Copy {
    /// Scans the sites where both sequences hold an unambiguous base.
    ///
    /// Returns `None` as soon as the mismatches seen so far exceed `max_mismatches`. Mismatches
    /// only grow as the scan continues, so this never discards a pair that could have passed.
    ///
    /// Implementations are kept out of line so every pair is a separate call. V8 cannot switch a
    /// running wasm function to optimized code, but later calls to a hot function use it, so a
    /// single long batch call still runs mostly optimized: the first `distances("tn93")` call
    /// over 1,500 x 1 kb took ~270 ms, about the same as later calls, against ~400 ms with
    /// baseline code only.
    fn count_planes(x: &[Word], y: &[Word], max_mismatches: f64) -> Option<Self>;

    /// Adds what one resolved ambiguous site contributes.
    fn add_resolved(&mut self, site: &CountsWithGc);

    /// Compared sites so far (weighted; an integer while every site weighs 1).
    fn sites(&self) -> f64;

    /// Adds the sites where a base faces an `N`, or two `N`s face each other (see
    /// [`FullWord`]): every site whose contribution is fixed by the base alone, counted with
    /// popcounts instead of one table lookup each. Returns how many sites it added.
    fn add_full_sites(&mut self, x: &Encoded, y: &Encoded, table: &Resolutions) -> u32;
}

/// Sites where one sequence holds a base and the other is compared as an `N`, by base, for each
/// side; and sites where both are compared as `N`.
struct FullSites {
    /// `[A, C, G, T]` held by sequence 1 facing an N in sequence 2.
    base1: [u32; 4],
    /// `[A, C, G, T]` held by sequence 2 facing an N in sequence 1.
    base2: [u32; 4],
    both: u32,
}

impl FullSites {
    /// Visits only the words where either sequence has such a site: ~5 popcounts per word for
    /// each side's. Only run for pairs where one sequence has any.
    #[inline(never)]
    fn count(x: &Encoded, y: &Encoded) -> FullSites {
        // [sites, lo, hi, lo & hi] per side: lo = C or T, hi = G or T, lo & hi = T.
        let (mut s1, mut s2) = ([0u32; 4], [0u32; 4]);
        let mut both = 0u32;
        // Bases of y facing x's Ns, and Ns facing Ns.
        for (w, [any1, n1]) in full_words(x, y) {
            let [lo2, hi2, v2] = y.words[w];
            let [any2, n2] = y.full.get(w).copied().unwrap_or_default();
            let m2 = v2 & any1;
            s2[0] += m2.count_ones();
            s2[1] += (m2 & lo2).count_ones();
            s2[2] += (m2 & hi2).count_ones();
            s2[3] += (m2 & lo2 & hi2).count_ones();
            both += (any1 & any2 & (n1 | n2)).count_ones();
        }
        // Bases of x facing y's Ns.
        for (w, [any2, _]) in full_words(y, x) {
            let [lo1, hi1, v1] = x.words[w];
            let m1 = v1 & any2;
            s1[0] += m1.count_ones();
            s1[1] += (m1 & lo1).count_ones();
            s1[2] += (m1 & hi1).count_ones();
            s1[3] += (m1 & lo1 & hi1).count_ones();
        }
        // A = n - (C + G + T) = n - (lo + hi - t); add t first so no step goes below zero.
        let by_base = |[n, lo, hi, t]: [u32; 4]| [n + t - lo - hi, lo - t, hi - t, t];
        FullSites { base1: by_base(s1), base2: by_base(s2), both }
    }

    /// Every site counted: a base on one side and an N on the other, or N on both.
    fn total(&self) -> u32 {
        self.base1.iter().sum::<u32>() + self.base2.iter().sum::<u32>() + self.both
    }
}

/// The non-zero words of `x`'s fully ambiguous plane that `other` also reaches, with their index.
#[inline]
fn full_words<'a>(x: &'a Encoded, other: &Encoded) -> impl Iterator<Item = (usize, FullWord)> + 'a {
    let len = other.words.len();
    x.full_index.iter().map(|&w| w as usize).take_while(move |&w| w < len).map(|w| (w, x.full[w]))
}

/// Adds every site of `full` to `tally`, scaled table entries per base. Returns the site count.
fn add_full_by_base(tally: &mut CountsWithGc, full: &FullSites, table: &Resolutions) -> u32 {
    for b in 0..4 {
        tally.add_scaled(table.get(1 << b, ANY), full.base1[b] as f64);
        tally.add_scaled(table.get(ANY, 1 << b), full.base2[b] as f64);
    }
    tally.add_scaled(table.get(ANY, ANY), full.both as f64);
    full.total()
}

impl Counts {
    /// The full bit-plane scan, ~10 popcounts per 64 sites (+1 for `gc1` with `GC`).
    #[inline(always)]
    fn scan_planes<const GC: bool>(x: &[Word], y: &[Word], max_mismatches: f64) -> Option<CountsWithGc> {
        let (mut sites, mut tv, mut ts, mut ts_y) = (0u32, 0u32, 0u32, 0u32);
        let (mut lo_n, mut hi_n, mut t_n, mut gc1) = (0u32, 0u32, 0u32, 0u32);
        let len = x.len().min(y.len());
        for (cx, cy) in x[..len].chunks(CHUNK).zip(y[..len].chunks(CHUNK)) {
            for (a, b) in cx.iter().zip(cy) {
                let [lo1, hi1, v1] = *a;
                let [lo2, hi2, v2] = *b;
                let v = v1 & v2;
                let tv_mask = (lo1 ^ lo2) & v; // purine <-> pyrimidine
                let ts_mask = (hi1 ^ hi2) & v & !tv_mask; // same class, different base
                sites += v.count_ones();
                tv += tv_mask.count_ones();
                ts += ts_mask.count_ones();
                ts_y += (ts_mask & lo1).count_ones();
                lo_n += (v & lo1).count_ones() + (v & lo2).count_ones();
                hi_n += (v & hi1).count_ones() + (v & hi2).count_ones();
                t_n += (v & lo1 & hi1).count_ones() + (v & lo2 & hi2).count_ones();
                if GC {
                    gc1 += (v & (lo1 ^ hi1)).count_ones(); // C = 01, G = 10
                }
            }
            if (ts + tv) as f64 > max_mismatches {
                return None;
            }
        }
        // lo = C or T, hi = G or T, lo & hi = T
        let (c, g, t) = (lo_n - t_n, hi_n - t_n, t_n);
        let a = 2 * sites - lo_n - g;
        let counts = Counts {
            sites: sites as f64,
            ts_purine: (ts - ts_y) as f64,
            ts_pyrimidine: ts_y as f64,
            transversions: tv as f64,
            freq: [a as f64, c as f64, g as f64, t as f64],
        };
        Some(CountsWithGc { counts, gc1: gc1 as f64 })
    }
}

impl Tally for Counts {
    #[inline(never)]
    fn count_planes(x: &[Word], y: &[Word], max_mismatches: f64) -> Option<Counts> {
        Counts::scan_planes::<false>(x, y, max_mismatches).map(|c| c.counts)
    }

    #[inline]
    fn add_resolved(&mut self, site: &CountsWithGc) {
        self.add_counts(&site.counts);
    }

    fn sites(&self) -> f64 {
        self.sites
    }

    fn add_full_sites(&mut self, x: &Encoded, y: &Encoded, table: &Resolutions) -> u32 {
        let mut c = CountsWithGc::default();
        let n = add_full_by_base(&mut c, &FullSites::count(x, y), table);
        self.add_counts(&c.counts);
        n
    }
}

impl Tally for CountsWithGc {
    #[inline(never)]
    fn count_planes(x: &[Word], y: &[Word], max_mismatches: f64) -> Option<CountsWithGc> {
        Counts::scan_planes::<true>(x, y, max_mismatches)
    }

    #[inline]
    fn add_resolved(&mut self, site: &CountsWithGc) {
        self.counts.add_counts(&site.counts);
        self.gc1 += site.gc1;
    }

    fn sites(&self) -> f64 {
        self.counts.sites
    }

    fn add_full_sites(&mut self, x: &Encoded, y: &Encoded, table: &Resolutions) -> u32 {
        add_full_by_base(self, &FullSites::count(x, y), table)
    }
}

impl Tally for Mismatches {
    /// 2 popcounts per 64 sites.
    #[inline(never)]
    fn count_planes(x: &[Word], y: &[Word], max_mismatches: f64) -> Option<Mismatches> {
        let (mut sites, mut mismatches) = (0u32, 0u32);
        let len = x.len().min(y.len());
        for (cx, cy) in x[..len].chunks(CHUNK).zip(y[..len].chunks(CHUNK)) {
            for (a, b) in cx.iter().zip(cy) {
                let [lo1, hi1, v1] = *a;
                let [lo2, hi2, v2] = *b;
                let v = v1 & v2;
                sites += v.count_ones();
                mismatches += (((lo1 ^ lo2) | (hi1 ^ hi2)) & v).count_ones();
            }
            if mismatches as f64 > max_mismatches {
                return None;
            }
        }
        Some(Mismatches { sites: sites as f64, mismatches: mismatches as f64 })
    }

    #[inline]
    fn add_resolved(&mut self, site: &CountsWithGc) {
        self.sites += site.counts.sites;
        self.mismatches += site.counts.mismatches();
    }

    /// Only needs how many sites have an N on one side or both: a base facing an N weighs the
    /// same whatever the base (checked in `Resolutions::new`).
    fn sites(&self) -> f64 {
        self.sites
    }

    fn add_full_sites(&mut self, x: &Encoded, y: &Encoded, table: &Resolutions) -> u32 {
        let (mut one, mut both) = (0u32, 0u32);
        for (w, [any1, n1]) in full_words(x, y) {
            let [any2, n2] = y.full.get(w).copied().unwrap_or_default();
            one += (y.words[w][2] & any1).count_ones();
            both += (any1 & any2 & (n1 | n2)).count_ones();
        }
        for (w, [any2, _]) in full_words(y, x) {
            one += (x.words[w][2] & any2).count_ones();
        }
        for (site, n) in [(table.get(1, ANY), one), (table.get(ANY, ANY), both)] {
            let n = n as f64;
            if n != 0.0 {
                self.sites += n * site.counts.sites;
                self.mismatches += n * site.counts.mismatches();
            }
        }
        one + both
    }
}

/// Sites the ambiguity passes added to a pair, as integers.
#[derive(Clone, Copy, Default)]
struct Visited {
    /// Compared sites (each weighs exactly 1).
    sites: u32,
    /// Sites tn93 counts toward `-g`'s ambiguity fraction in resolve mode (`ambig_count` in
    /// `computeTN93`): a base facing a code that contains it, or two ambiguity codes, whether or
    /// not they share a base. A base facing a code without it is averaged and does not count.
    ambiguous: u32,
}

/// Whether tn93 counts a site holding codes `c1` and `c2` (at least one ambiguous) as resolved.
#[inline]
fn counts_as_ambiguous(c1: u8, c2: u8) -> bool {
    c1 & c2 != 0 || (c1.count_ones() > 1 && c2.count_ones() > 1)
}

/// Ambiguity fallback: adds the sites where at least one sequence holds a partial ambiguity
/// code and the other holds any base, a partial code or (in modes that compare them) an `N`.
/// Only called when one of the two tracks ambiguities. Visits only the (sparse) ambiguous sites
/// of each sequence. Counts what it visits only with `COUNT` (for a maximum ambiguity fraction):
/// always counting cost the other modes ~8%.
fn add_ambiguous_sites<T: Tally, const COUNT: bool>(x: &Encoded, y: &Encoded, table: &Resolutions, tally: &mut T) -> Visited {
    let mut visited = Visited::default();
    let mut add = |c1: u8, c2: u8| {
        if c1 != 0 && c2 != 0 {
            tally.add_resolved(table.get(c1, c2));
            if COUNT {
                visited.sites += 1;
                visited.ambiguous += counts_as_ambiguous(c1, c2) as u32;
            }
        }
    };
    if let Some(ax) = &x.ambiguity {
        for (&k, &cx) in ax.positions.iter().zip(&ax.codes) {
            add(cx, y.code_at(k as usize));
        }
    }
    if let Some(ay) = &y.ambiguity {
        for (&k, &cy) in ay.positions.iter().zip(&ay.codes) {
            // Sites ambiguous in both were counted above; base_code_at gives 0 for them.
            add(x.base_code_at(k as usize), cy);
        }
    }
    visited
}

/// How a set compares ambiguity codes: the resolution table for its mode and, for tn93's
/// `-a resolve -g FRACTION`, the averaging table it falls back to for too ambiguous pairs.
pub(crate) struct Comparison {
    pub table: Resolutions,
    /// `(fraction, average)`: a pair whose resolved ambiguous sites (see [`Visited`]) are at
    /// least `fraction` of its compared sites is averaged instead. Only with `fraction < 1`.
    pub average_above: Option<(f64, Resolutions)>,
}

impl Comparison {
    pub fn new(mode: Ambiguity, max_ambiguity_fraction: Option<f64>) -> Self {
        let average_above = max_ambiguity_fraction
            .filter(|&f| f < 1.0)
            .map(|f| (f, Resolutions::new(Ambiguity::Average)));
        Comparison { table: Resolutions::new(mode), average_above }
    }
}

/// Scans one pair of encoded sequences: fast kernel, plus the passes over sites compared as `N`
/// and over partial ambiguity codes when either sequence has any (which only happens when
/// ambiguity codes are compared at all).
///
/// Early exit only looks at the fast kernel's mismatches. The other passes only add mismatches,
/// so a pair that already has too many cannot come back under the bound, whichever table the
/// pair ends up using: resolving and averaging compare the same sites.
#[inline]
pub(crate) fn count_pair<T: Tally>(x: &Encoded, y: &Encoded, comparison: Option<&Comparison>, max_mismatches: f64) -> Option<T> {
    let mut tally = T::count_planes(&x.words, &y.words, max_mismatches)?;
    if let Some(comparison) = comparison {
        if x.has_ambiguity() || y.has_ambiguity() {
            add_ambiguity_passes(x, y, comparison, &mut tally);
        }
    }
    Some(tally)
}

/// The passes over compared ambiguity codes, in one out-of-line call, so that `count_pair`,
/// which is inlined into every row loop, stays small. With a maximum ambiguity fraction, the
/// passes run with the resolution table and, if the pair turns out too ambiguous, again with the
/// averaging table on top of the same fast-kernel counts.
#[inline(never)]
fn add_ambiguity_passes<T: Tally>(x: &Encoded, y: &Encoded, comparison: &Comparison, tally: &mut T) {
    let Some((fraction, average)) = &comparison.average_above else {
        add_passes::<T, false>(x, y, &comparison.table, tally);
        return;
    };
    let planes = *tally;
    let visited = add_passes::<T, true>(x, y, &comparison.table, tally);
    // Exact integer site counts: the weighted sum can round below an integer.
    let sites = planes.sites() + visited.sites as f64;
    if sites * fraction <= visited.ambiguous as f64 {
        *tally = planes;
        add_passes::<T, false>(x, y, average, tally);
    }
}

fn add_passes<T: Tally, const COUNT: bool>(x: &Encoded, y: &Encoded, table: &Resolutions, tally: &mut T) -> Visited {
    let mut visited = Visited::default();
    if x.has_full() || y.has_full() {
        // In the modes that average too ambiguous pairs (resolve), these are bases facing an N,
        // or two Ns: every one counts as ambiguous.
        let n = tally.add_full_sites(x, y, table);
        visited = Visited { sites: n, ambiguous: n };
    }
    if x.ambiguity.is_some() || y.ambiguity.is_some() {
        let sparse = add_ambiguous_sites::<T, COUNT>(x, y, table, tally);
        visited.sites += sparse.sites;
        visited.ambiguous += sparse.ambiguous;
    }
    visited
}
