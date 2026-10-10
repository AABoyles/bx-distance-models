//! Sequence encodings.
//!
//! Every sequence gets the fast encoding: two bit-planes plus a validity mask, 64 sites per
//! `u64` word (A=00 C=01 G=10 T=11, so the low bit marks pyrimidines). Sites that are not an
//! unambiguous base have their validity bit cleared and are skipped by the fast kernel.
//!
//! When ambiguity resolution is requested, sequences that contain partial IUPAC ambiguity codes
//! (R, Y, S, W, K, M, B, D, H, V) additionally keep a sparse list of those sites with their 4-bit
//! IUPAC codes (NtSeq-style one-bit-per-base sets) for the slower fallback kernel.

/// IUPAC nucleotide as a 4-bit set: A=1, C=2, G=4, T=8 (`N` = all four).
/// 0 = not a nucleotide (gap, `?`, anything unrecognized). Case-insensitive; U is read as T.
pub(crate) static IUPAC: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes: [(u8, u8); 16] = [
        (b'A', 1), (b'C', 2), (b'G', 4), (b'T', 8), (b'U', 8),
        (b'R', 1 | 4), (b'Y', 2 | 8), (b'S', 2 | 4), (b'W', 1 | 8), (b'K', 4 | 8), (b'M', 1 | 2),
        (b'B', 2 | 4 | 8), (b'D', 1 | 4 | 8), (b'H', 1 | 2 | 8), (b'V', 1 | 2 | 4), (b'N', 15),
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

/// How sites holding partial ambiguity codes (R, Y, ...) are compared.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ambiguity {
    /// Skip any site where either sequence is ambiguous (fast path only).
    Skip,
    /// Count a site as a match if the two codes share a base; otherwise average.
    Resolve,
    /// Average over every pair of bases the two codes allow.
    Average,
}

impl Ambiguity {
    pub fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "skip" => Some(Ambiguity::Skip),
            "resolve" => Some(Ambiguity::Resolve),
            "average" => Some(Ambiguity::Average),
            _ => None,
        }
    }
}

/// 64 sites of the fast encoding: `[lo, hi, valid]`.
pub(crate) type Word = [u64; 3];

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
    /// Sites that can take part in a comparison: unambiguous bases, plus partial ambiguity codes
    /// when they are tracked. Upper bound on the compared sites of any pair with this sequence.
    pub informative: u32,
    /// Present only when ambiguities are tracked and the sequence contains at least one.
    pub ambiguity: Option<AmbiguityCodes>,
}

impl Encoded {
    pub fn new(seq: &[u8], track_ambiguity: bool) -> Self {
        let mut words = Vec::with_capacity(seq.len().div_ceil(64));
        let mut mask = Vec::new();
        let mut informative = 0u32;
        let mut any_ambiguous = false;
        for chunk in seq.chunks(64) {
            let (mut lo, mut hi, mut valid, mut amb) = (0u64, 0u64, 0u64, 0u64);
            for (k, &ch) in chunk.iter().enumerate() {
                let base = BASE[ch as usize] as u64;
                if base < 4 {
                    valid |= 1 << k;
                    lo |= (base & 1) << k;
                    hi |= (base >> 1) << k;
                } else if !matches!(IUPAC[ch as usize], 0 | 15) {
                    amb |= 1 << k;
                }
            }
            words.push([lo, hi, valid]);
            if track_ambiguity {
                mask.push(amb);
                any_ambiguous |= amb != 0;
            }
            informative += valid.count_ones();
        }

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

        Encoded { words, informative, ambiguity }
    }

    /// IUPAC code of site `k` if it holds an unambiguous base, else 0.
    #[inline]
    pub fn base_code_at(&self, k: usize) -> u8 {
        let Some(&[lo, hi, valid]) = self.words.get(k >> 6) else { return 0 };
        let bit = k & 63;
        if (valid >> bit) & 1 == 0 {
            return 0;
        }
        1 << (((lo >> bit) & 1) | (((hi >> bit) & 1) << 1))
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
