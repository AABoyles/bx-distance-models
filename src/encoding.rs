//! Sequence encodings.
//!
//! Every sequence gets the fast encoding: two bit-planes plus a validity mask, 64 sites per
//! `u64` word (A=00 C=01 G=10 T=11, so the low bit marks pyrimidines). Sites that are not an
//! unambiguous base have their validity bit cleared and are skipped by the fast kernel.
//!
//! When ambiguity resolution is requested, sequences that contain partial IUPAC ambiguity codes
//! (R, Y, S, W, K, M, B, D, H, V) additionally keep a sparse list of those sites with their 4-bit
//! IUPAC codes (NtSeq-style one-bit-per-base sets) for the slower fallback kernel.
//!
//! Sites that stand for any base (`N` and `?` where the mode compares them, and gaps in GapMM
//! mode) are too common in real alignments for the sparse list. They get their own bit-plane
//! instead: a base facing one contributes a fixed quarter to each of its four (base, X) cells, so
//! a second popcount pass counts them without visiting each site.

/// IUPAC nucleotide as a 4-bit set: A=1, C=2, G=4, T=8 (`N` and `?` = all four, as in tn93).
/// 0 = not a nucleotide (gap, `.`, anything unrecognized). Case-insensitive; U is read as T.
pub(crate) static IUPAC: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes: [(u8, u8); 17] = [
        (b'A', 1), (b'C', 2), (b'G', 4), (b'T', 8), (b'U', 8),
        (b'R', 1 | 4), (b'Y', 2 | 8), (b'S', 2 | 4), (b'W', 1 | 8), (b'K', 4 | 8), (b'M', 1 | 2),
        (b'B', 2 | 4 | 8), (b'D', 1 | 4 | 8), (b'H', 1 | 2 | 8), (b'V', 1 | 2 | 4), (b'N', 15), (b'?', 15),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (ch, code) = codes[i];
        t[ch as usize] = code;
        t[ch.to_ascii_lowercase() as usize] = code;
        i += 1;
    }
    t
};

/// 2-bit base index (A=0 C=1 G=2 T=3) of an unambiguous base, 4 for anything else.
/// The base index is also the bit position of the base in its IUPAC code.
pub(crate) static BASE: [u8; 256] = {
    let mut t = [4u8; 256];
    let mut b = 0;
    while b < 256 {
        t[b] = match IUPAC[b] { 1 => 0, 2 => 1, 4 => 2, 8 => 3, _ => 4 };
        b += 1;
    }
    t
};

/// How sites holding ambiguity codes are compared; the modes of HIV-TRACE's `tn93 -a`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ambiguity {
    /// Skip any site where either sequence is ambiguous (fast path only).
    Skip,
    /// Count a site as a match if the two codes share a base; otherwise average.
    Resolve,
    /// Average over every pair of bases the two codes allow.
    Average,
    /// As `Average`, and a gap facing anything but a gap is compared as an `N`.
    GapMm,
}

impl Ambiguity {
    pub fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "skip" => Some(Ambiguity::Skip),
            "resolve" => Some(Ambiguity::Resolve),
            "average" => Some(Ambiguity::Average),
            "gapmm" => Some(Ambiguity::GapMm),
            _ => None,
        }
    }

    /// Whether sites holding ambiguity codes are compared at all.
    pub fn tracks_ambiguity(self) -> bool {
        self != Ambiguity::Skip
    }

    /// Whether `N` and `?` are compared as fourfold ambiguity codes: in every mode that
    /// compares ambiguity codes, as tn93 does.
    pub fn compares_n(self) -> bool {
        self != Ambiguity::Skip
    }

    /// Whether a gap (or any other character that is not a nucleotide) facing a non-gap is
    /// compared as an `N`.
    pub fn gap_is_n(self) -> bool {
        self == Ambiguity::GapMm
    }
}

/// IUPAC code standing for any base.
pub(crate) const ANY: u8 = 15;

/// 64 sites of the fast encoding: `[lo, hi, valid]`.
pub(crate) type Word = [u64; 3];

/// 64 sites of the fully ambiguous plane: `[any, n]`. `any` marks every site compared as an `N`
/// (`N`, `?`, and in GapMM mode gaps and other non-nucleotides); `n` only actual `N` and `?`.
/// Two gaps never face each other as `N`s, so a pair compares `any1 & any2 & (n1 | n2)`.
pub(crate) type FullWord = [u64; 2];

/// The partial-ambiguity sites of a sequence that contains any.
pub(crate) struct AmbiguityCodes {
    /// Sites holding a 2- or 3-base ambiguity code, 64 per word (for O(1) lookups).
    pub mask: Vec<u64>,
    /// The same sites in ascending order, with their 4-bit IUPAC codes.
    pub positions: Vec<u32>,
    pub codes: Vec<u8>,
}

pub(crate) struct Encoded {
    pub words: Vec<Word>,
    /// Sites that can take part in a comparison: unambiguous bases, plus every site the mode
    /// compares as an ambiguity code. Upper bound on the compared sites of any pair with this
    /// sequence.
    pub informative: u32,
    /// Present only when ambiguities are tracked and the sequence contains at least one.
    pub ambiguity: Option<AmbiguityCodes>,
    /// The fully ambiguous plane, one entry per word; empty unless the mode compares such sites
    /// and the sequence has at least one.
    pub full: Vec<FullWord>,
    /// Indices of the words of `full` that are not zero, ascending. Gaps and N cluster (ragged
    /// ends, runs), so a pair only visits these.
    pub full_index: Vec<u32>,
}

impl Encoded {
    pub fn new(seq: &[u8], mode: Ambiguity) -> Self {
        let (track, compare_n, gap_is_n) = (mode.tracks_ambiguity(), mode.compares_n(), mode.gap_is_n());
        let mut words = Vec::with_capacity(seq.len().div_ceil(64));
        let mut mask = Vec::new();
        let mut full = Vec::new();
        let mut informative = 0u32;
        let (mut any_ambiguous, mut any_full) = (false, false);
        for chunk in seq.chunks(64) {
            let (mut lo, mut hi, mut valid, mut amb, mut any, mut n) = (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
            for (k, &ch) in chunk.iter().enumerate() {
                let base = BASE[ch as usize] as u64;
                if base < 4 {
                    valid |= 1 << k;
                    lo |= (base & 1) << k;
                    hi |= (base >> 1) << k;
                } else {
                    match IUPAC[ch as usize] {
                        ANY => n |= (compare_n as u64) << k,
                        0 => any |= (gap_is_n as u64) << k,
                        _ => amb |= 1 << k,
                    }
                }
            }
            any |= n;
            words.push([lo, hi, valid]);
            if track {
                mask.push(amb);
                any_ambiguous |= amb != 0;
                full.push([any, n]);
                any_full |= any != 0;
            }
            informative += valid.count_ones() + any.count_ones();
        }
        if !any_full {
            full = Vec::new();
        }
        let full_index = (0..full.len() as u32).filter(|&w| full[w as usize][0] != 0).collect();

        let ambiguity = if any_ambiguous {
            let (mut positions, mut codes) = (Vec::new(), Vec::new());
            for (w, &m) in mask.iter().enumerate() {
                let mut m = m;
                while m != 0 {
                    let k = w * 64 + m.trailing_zeros() as usize;
                    m &= m - 1;
                    positions.push(k as u32);
                    codes.push(IUPAC[seq[k] as usize]);
                }
            }
            informative += positions.len() as u32;
            Some(AmbiguityCodes { mask, positions, codes })
        } else {
            None
        };

        Encoded { words, informative, ambiguity, full, full_index }
    }

    /// IUPAC code of site `k` if it holds an unambiguous base, [`ANY`] if the mode compares it
    /// as an `N`, else 0. Partial ambiguity codes give 0 too.
    #[inline]
    pub fn base_code_at(&self, k: usize) -> u8 {
        let Some(&[lo, hi, valid]) = self.words.get(k >> 6) else { return 0 };
        let bit = k & 63;
        if (valid >> bit) & 1 == 0 {
            let any = self.full.get(k >> 6).is_some_and(|f| (f[0] >> bit) & 1 == 1);
            return if any { ANY } else { 0 };
        }
        1 << (((lo >> bit) & 1) | (((hi >> bit) & 1) << 1))
    }

    /// Whether any site is compared as an `N` (the pass over the fully ambiguous plane is needed).
    #[inline]
    pub fn has_full(&self) -> bool {
        !self.full.is_empty()
    }

    /// Whether any site holds an ambiguity code the mode compares (pairs with this sequence
    /// need more than the fast kernel).
    #[inline]
    pub fn has_ambiguity(&self) -> bool {
        self.ambiguity.is_some() || self.has_full()
    }

    /// Whether site `k` holds a tracked partial ambiguity code.
    #[inline]
    pub fn is_ambiguous_at(&self, k: usize) -> bool {
        self.ambiguity.as_ref().is_some_and(|a| a.mask.get(k >> 6).is_some_and(|m| (m >> (k & 63)) & 1 == 1))
    }

    /// IUPAC code of site `k` (a base or a tracked ambiguity code), or 0 if it is not compared.
    #[inline]
    pub fn code_at(&self, k: usize) -> u8 {
        if self.is_ambiguous_at(k) {
            let a = self.ambiguity.as_ref().unwrap();
            return a.codes[a.positions.binary_search(&(k as u32)).unwrap()];
        }
        self.base_code_at(k)
    }
}
