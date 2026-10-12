//! Per-pair site counts and the kernels that produce them.
//!
//! Every metric in this crate is a function of [`Counts`], so each pair is scanned once no
//! matter how many metrics are wanted. Metrics that only need compared sites and mismatches
//! (SNPs, Hamming, p-distance, JC69) can use the cheaper [`Mismatches`] scan instead.
use crate::encoding::{Ambiguity, Encoded, Word, BASE};

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
        self.sites += other.sites;
        self.ts_purine += other.ts_purine;
        self.ts_pyrimidine += other.ts_pyrimidine;
        self.transversions += other.transversions;
        for (f, o) in self.freq.iter_mut().zip(other.freq) {
            *f += o;
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
        // Codes 0 (not a base) and 15 (N) are never compared.
        for c1 in 1..15u8 {
            for c2 in 1..15u8 {
                table[(c1 as usize) << 4 | c2 as usize].add_codes(c1, c2, 1.0, mode);
            }
        }
        Resolutions(table)
    }

    #[inline]
    fn get(&self, c1: u8, c2: u8) -> &CountsWithGc {
        &self.0[(c1 as usize) << 4 | c2 as usize]
    }
}

/// Mismatches are checked for early exit every 8 words (512 sites).
const CHUNK: usize = 8;

/// What a scan of one pair accumulates: full [`Counts`], or only [`Mismatches`].
pub(crate) trait Tally: Sized {
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
}

/// Ambiguity fallback: adds the sites where at least one sequence holds a partial ambiguity
/// code and the other holds any base. Only called when one of the two tracks ambiguities.
/// Visits only the (sparse) ambiguous sites of each sequence.
pub(crate) fn add_ambiguous_sites<T: Tally>(x: &Encoded, y: &Encoded, table: &Resolutions, tally: &mut T) {
    if let Some(ax) = &x.ambiguity {
        for (&k, &cx) in ax.positions.iter().zip(&ax.codes) {
            let cy = y.code_at(k as usize);
            if cy != 0 {
                tally.add_resolved(table.get(cx, cy));
            }
        }
    }
    if let Some(ay) = &y.ambiguity {
        for (&k, &cy) in ay.positions.iter().zip(&ay.codes) {
            // Sites ambiguous in both were counted above.
            let cx = x.base_code_at(k as usize);
            if cx != 0 {
                tally.add_resolved(table.get(cx, cy));
            }
        }
    }
}

/// Scans one pair of encoded sequences: fast kernel, plus the ambiguity fallback when either
/// sequence carries ambiguity codes (which only happens when `table` is set).
#[inline]
pub(crate) fn count_pair<T: Tally>(x: &Encoded, y: &Encoded, table: Option<&Resolutions>, max_mismatches: f64) -> Option<T> {
    let mut tally = T::count_planes(&x.words, &y.words, max_mismatches)?;
    if let Some(table) = table {
        if x.ambiguity.is_some() || y.ambiguity.is_some() {
            add_ambiguous_sites(x, y, table, &mut tally);
        }
    }
    Some(tally)
}
