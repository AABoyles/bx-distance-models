//! Kernel benchmark: the repo's current byte-loop kernels vs. pre-encoded alternatives.
//! Raw `extern "C"` exports (no wasm-bindgen) so we time the kernels, not the glue.
use std::alloc::{alloc as sys_alloc, Layout};
use std::slice;

#[no_mangle]
pub extern "C" fn alloc(size: usize) -> *mut u8 {
    unsafe { sys_alloc(Layout::from_size_align(size.max(16), 16).unwrap()) }
}

// =====================================================================
// 1. CURRENT: verbatim logic from src/metrics/{basic,evolutionary}.rs
// =====================================================================
#[inline] fn is_valid_base(b: u8) -> bool { matches!(b, b'A' | b'C' | b'G' | b'T') }
#[inline] fn is_purine(b: u8) -> bool { matches!(b, b'A' | b'G') }
#[inline] fn is_pyrimidine(b: u8) -> bool { matches!(b, b'C' | b'T') }

pub fn snps_current(bytes1: &[u8], bytes2: &[u8]) -> f64 {
    let len = bytes1.len().min(bytes2.len());
    if len == 0 { return 0.0; }
    let mut mismatches = 0u32;
    for i in 0..len {
        let a = bytes1[i]; let b = bytes2[i];
        if a == b'-' || b == b'-' || a == b'N' || b == b'N' { continue; }
        if a != b { mismatches += 1; }
    }
    mismatches as f64
}

pub fn hamming_current(bytes1: &[u8], bytes2: &[u8]) -> f64 {
    let len = bytes1.len().min(bytes2.len());
    if len == 0 { return 1.0; }
    let mut mismatches = 0u32;
    let mut compared = 0u32;
    for i in 0..len {
        let a = bytes1[i]; let b = bytes2[i];
        if a == b'-' || b == b'-' || a == b'N' || b == b'N' { continue; }
        compared += 1;
        if a != b { mismatches += 1; }
    }
    if compared == 0 { 1.0 } else { mismatches as f64 / compared as f64 }
}

pub fn tn93_current(bytes1: &[u8], bytes2: &[u8]) -> f64 {
    let len = bytes1.len().min(bytes2.len());
    if len == 0 { return 1.0; }
    let mut total_sites = 0u32;
    let mut transitions1 = 0u32;
    let mut transitions2 = 0u32;
    let mut transversions = 0u32;
    let mut freq_a = 0u32; let mut freq_c = 0u32; let mut freq_g = 0u32; let mut freq_t = 0u32;
    for i in 0..len {
        let a = bytes1[i]; let b = bytes2[i];
        if !is_valid_base(a) || !is_valid_base(b) { continue; }
        total_sites += 1;
        match a { b'A' => freq_a += 1, b'C' => freq_c += 1, b'G' => freq_g += 1, b'T' => freq_t += 1, _ => {} }
        match b { b'A' => freq_a += 1, b'C' => freq_c += 1, b'G' => freq_g += 1, b'T' => freq_t += 1, _ => {} }
        if a == b { continue; }
        if is_purine(a) && is_purine(b) { transitions1 += 1; }
        else if is_pyrimidine(a) && is_pyrimidine(b) { transitions2 += 1; }
        else { transversions += 1; }
    }
    if total_sites == 0 { return 1.0; }
    let total2 = (total_sites * 2) as f64;
    let g_a = freq_a as f64 / total2; let g_c = freq_c as f64 / total2;
    let g_g = freq_g as f64 / total2; let g_t = freq_t as f64 / total2;
    let g_r = g_a + g_g; let g_y = g_c + g_t;
    if g_r == 0.0 || g_y == 0.0 { return hamming_current(bytes1, bytes2); }
    let p1 = transitions1 as f64 / total_sites as f64;
    let p2 = transitions2 as f64 / total_sites as f64;
    let q = transversions as f64 / total_sites as f64;
    let a1 = (2.0 * g_a * g_g) / g_r;
    let a2 = (2.0 * g_c * g_t) / g_y;
    let b = 2.0 * (g_r * g_y - (g_a * g_g * g_y) / g_r - (g_c * g_t * g_r) / g_y);
    let term1 = 1.0 - p1 / a1 - q / (2.0 * g_r);
    let term2 = 1.0 - p2 / a2 - q / (2.0 * g_y);
    let term3 = 1.0 - q / (2.0 * g_r * g_y);
    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 { return hamming_current(bytes1, bytes2); }
    let d = -a1 * term1.ln() - a2 * term2.ln() - b * term3.ln();
    d.max(0.0).min(1.0)
}

// Same TN93 formula, fed from counts. The hamming fallback is mism/sites, which equals the
// repo's hamming_distance() for ACGTN- input (the benchmark alphabet).
#[inline]
fn tn93_finish(sites: u32, ts_r: u32, ts_y: u32, tv: u32, fa: u32, fc: u32, fg: u32, ft: u32) -> f64 {
    if sites == 0 { return 1.0; }
    let total2 = (sites * 2) as f64;
    let g_a = fa as f64 / total2; let g_c = fc as f64 / total2;
    let g_g = fg as f64 / total2; let g_t = ft as f64 / total2;
    let g_r = g_a + g_g; let g_y = g_c + g_t;
    let ham = (ts_r + ts_y + tv) as f64 / sites as f64;
    if g_r == 0.0 || g_y == 0.0 { return ham; }
    let p1 = ts_r as f64 / sites as f64;
    let p2 = ts_y as f64 / sites as f64;
    let q = tv as f64 / sites as f64;
    let a1 = (2.0 * g_a * g_g) / g_r;
    let a2 = (2.0 * g_c * g_t) / g_y;
    let b = 2.0 * (g_r * g_y - (g_a * g_g * g_y) / g_r - (g_c * g_t * g_r) / g_y);
    let term1 = 1.0 - p1 / a1 - q / (2.0 * g_r);
    let term2 = 1.0 - p2 / a2 - q / (2.0 * g_y);
    let term3 = 1.0 - q / (2.0 * g_r * g_y);
    if term1 <= 0.0 || term2 <= 0.0 || term3 <= 0.0 { return ham; }
    let d = -a1 * term1.ln() - a2 * term2.ln() - b * term3.ln();
    d.max(0.0).min(1.0)
}

// =====================================================================
// 2. LUT: pre-encode once to 1 byte/site (A0 C1 G2 T3 other4), then a
//    branch-free 5x5 count matrix per pair (HIV-TRACE tn93 style).
// =====================================================================
static CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0; t[b'C' as usize] = 1; t[b'G' as usize] = 2; t[b'T' as usize] = 3;
    t
};

fn tn93_lut(c1: &[u8], c2: &[u8]) -> f64 {
    let mut m = [0u32; 32];
    for (&a, &b) in c1.iter().zip(c2) {
        m[(a as usize * 5 + b as usize) & 31] += 1;
    }
    let at = |i: usize, j: usize| m[i * 5 + j];
    let (mut sites, mut diag) = (0u32, 0u32);
    let mut f = [0u32; 4];
    for i in 0..4 {
        for j in 0..4 {
            let c = at(i, j);
            sites += c; f[i] += c; f[j] += c;
            if i == j { diag += c; }
        }
    }
    let ts_r = at(0, 2) + at(2, 0);
    let ts_y = at(1, 3) + at(3, 1);
    let tv = sites - diag - ts_r - ts_y;
    tn93_finish(sites, ts_r, ts_y, tv, f[0], f[1], f[2], f[3])
}

// =====================================================================
// 3. BITS: 2 bit-planes + validity mask, 64 sites per u64 word.
//    A=00 C=01 G=10 T=11 (lo bit = pyrimidine). Layout: words x [lo, hi, valid].
//    transversion = lo differs; transition = hi differs & lo same.
// =====================================================================
type W = [u64; 3];

#[derive(Clone, Copy, Default)]
#[repr(C)]
pub struct Tot { v: u32, y: u32, h: u32, t: u32 } // valid, C+T, G+T, T (over own valid sites)

fn encode_bits(src: &[u8], dst: &mut [W]) -> Tot {
    let mut tot = Tot::default();
    for (w, chunk) in src.chunks(64).enumerate() {
        let (mut lo, mut hi, mut v) = (0u64, 0u64, 0u64);
        for (k, &b) in chunk.iter().enumerate() {
            let c = CODE[b as usize] as u64;
            if c < 4 { v |= 1 << k; lo |= (c & 1) << k; hi |= (c >> 1) << k; }
        }
        dst[w] = [lo, hi, v];
        tot.v += v.count_ones(); tot.y += lo.count_ones(); tot.h += hi.count_ones(); tot.t += (lo & hi).count_ones();
    }
    tot
}

/// Straightforward: 10 popcounts per 64 sites, no precomputation.
fn tn93_bits(x: &[W], y: &[W]) -> f64 {
    let (mut sites, mut tsr, mut tsy, mut tv, mut ny, mut nh, mut nt) = (0u32, 0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    for (a, b) in x.iter().zip(y) {
        let [lo1, hi1, v1] = *a; let [lo2, hi2, v2] = *b;
        let v = v1 & v2;
        let dlo = (lo1 ^ lo2) & v;
        let ts = (hi1 ^ hi2) & v & !dlo;
        sites += v.count_ones(); tv += dlo.count_ones();
        tsy += (ts & lo1).count_ones(); tsr += (ts & !lo1).count_ones();
        ny += (v & lo1).count_ones() + (v & lo2).count_ones();
        nh += (v & hi1).count_ones() + (v & hi2).count_ones();
        nt += (v & lo1 & hi1).count_ones() + (v & lo2 & hi2).count_ones();
    }
    let (fc, fg, ft) = (ny - nt, nh - nt, nt);
    let fa = 2 * sites - ny - fg;
    tn93_finish(sites, tsr, tsy, tv, fa, fc, fg, ft)
}

/// Faster: per-sequence base totals precomputed once; per pair only subtract the
/// (rare) sites valid in exactly one sequence. 4 popcounts per 64 sites in the common case.
#[inline]
fn tn93_bits_fast_inner(x: &[W], tx: Tot, y: &[W], ty: Tot, thr: f64) -> f64 {
    let (mut sites, mut ts, mut tsy, mut tv) = (0u32, 0u32, 0u32, 0u32);
    let (mut ry, mut rh, mut rt) = (0u32, 0u32, 0u32);
    let bound = thr * tx.v.min(ty.v) as f64; // compared sites <= min(valid1, valid2)
    for (cx, cy) in x.chunks(8).zip(y.chunks(8)) {
        for (a, b) in cx.iter().zip(cy) {
            let [lo1, hi1, v1] = *a; let [lo2, hi2, v2] = *b;
            let v = v1 & v2;
            let dlo = (lo1 ^ lo2) & v;
            let t = (hi1 ^ hi2) & v & !dlo;
            sites += v.count_ones(); tv += dlo.count_ones();
            ts += t.count_ones(); tsy += (t & lo1).count_ones();
            let only = (v1 | v2) & !v;
            if only != 0 {
                let (m1, m2) = (only & v1, only & v2);
                ry += (m1 & lo1).count_ones() + (m2 & lo2).count_ones();
                rh += (m1 & hi1).count_ones() + (m2 & hi2).count_ones();
                rt += (m1 & lo1 & hi1).count_ones() + (m2 & lo2 & hi2).count_ones();
            }
        }
        // p_final >= mism_so_far / min(V1,V2), and TN93 >= p, so this pair can never pass.
        if (ts + tv) as f64 > bound { return f64::INFINITY; }
    }
    let ny = tx.y + ty.y - ry; let nh = tx.h + ty.h - rh; let nt = tx.t + ty.t - rt;
    let (fc, fg, ft) = (ny - nt, nh - nt, nt);
    let fa = 2 * sites - ny - fg;
    tn93_finish(sites, ts - tsy, tsy, tv, fa, fc, fg, ft)
}

fn hamming_bits(x: &[W], y: &[W]) -> (u32, u32) {
    let (mut sites, mut mism) = (0u32, 0u32);
    for (a, b) in x.iter().zip(y) {
        let [lo1, hi1, v1] = *a; let [lo2, hi2, v2] = *b;
        let v = v1 & v2;
        sites += v.count_ones();
        mism += (((lo1 ^ lo2) | (hi1 ^ hi2)) & v).count_ones();
    }
    (sites, mism)
}

// =====================================================================
// Pairwise drivers (upper triangle, row-major, same order as the repo)
// =====================================================================
unsafe fn ascii<'a>(p: *const u8, n: usize, l: usize) -> Vec<&'a [u8]> {
    (0..n).map(|i| slice::from_raw_parts(p.add(i * l), l)).collect()
}
unsafe fn bits<'a>(p: *const W, n: usize, words: usize) -> Vec<&'a [W]> {
    (0..n).map(|i| slice::from_raw_parts(p.add(i * words), words)).collect()
}

#[no_mangle]
pub unsafe extern "C" fn pairs_tn93_current(p: *const u8, n: usize, l: usize, out: *mut f64) {
    let s = ascii(p, n, l); let mut k = 0;
    for i in 0..n { for j in i + 1..n { *out.add(k) = tn93_current(s[i], s[j]); k += 1; } }
}
#[no_mangle]
pub unsafe extern "C" fn pairs_hamming_current(p: *const u8, n: usize, l: usize, out: *mut f64) {
    let s = ascii(p, n, l); let mut k = 0;
    for i in 0..n { for j in i + 1..n { *out.add(k) = hamming_current(s[i], s[j]); k += 1; } }
}
/// Mirrors compute_pairwise_distances_flat_all: three separate passes per pair.
#[no_mangle]
pub unsafe extern "C" fn pairs_all_current(p: *const u8, n: usize, l: usize, out: *mut f64) {
    let s = ascii(p, n, l); let np = n * (n - 1) / 2; let mut k = 0;
    for i in 0..n { for j in i + 1..n {
        *out.add(k) = hamming_current(s[i], s[j]);
        *out.add(k + np) = tn93_current(s[i], s[j]);
        *out.add(k + 2 * np) = snps_current(s[i], s[j]);
        k += 1;
    } }
}

#[no_mangle]
pub unsafe extern "C" fn encode_lut(p: *const u8, len: usize, out: *mut u8) {
    let src = slice::from_raw_parts(p, len); let dst = slice::from_raw_parts_mut(out, len);
    for (d, &s) in dst.iter_mut().zip(src) { *d = CODE[s as usize]; }
}
#[no_mangle]
pub unsafe extern "C" fn pairs_tn93_lut(p: *const u8, n: usize, l: usize, out: *mut f64) {
    let s = ascii(p, n, l); let mut k = 0;
    for i in 0..n { for j in i + 1..n { *out.add(k) = tn93_lut(s[i], s[j]); k += 1; } }
}

#[no_mangle]
pub unsafe extern "C" fn encode_bits_all(p: *const u8, n: usize, l: usize, out: *mut W, tots: *mut Tot) {
    let words = (l + 63) / 64;
    for i in 0..n {
        let src = slice::from_raw_parts(p.add(i * l), l);
        let dst = slice::from_raw_parts_mut(out.add(i * words), words);
        *tots.add(i) = encode_bits(src, dst);
    }
}
#[no_mangle]
pub unsafe extern "C" fn pairs_tn93_bits(p: *const W, n: usize, words: usize, out: *mut f64) {
    let s = bits(p, n, words); let mut k = 0;
    for i in 0..n { for j in i + 1..n { *out.add(k) = tn93_bits(s[i], s[j]); k += 1; } }
}
#[no_mangle]
pub unsafe extern "C" fn pairs_tn93_bits_fast(p: *const W, t: *const Tot, n: usize, words: usize, out: *mut f64) {
    pairs_tn93_bits_rows(p, t, n, words, 0, n, f64::INFINITY, out);
}
/// Rows [i0, i1) only — used to split the triangle across workers. Returns #pairs <= thr.
#[no_mangle]
pub unsafe extern "C" fn pairs_tn93_bits_rows(p: *const W, t: *const Tot, n: usize, words: usize,
                                              i0: usize, i1: usize, thr: f64, out: *mut f64) -> u32 {
    let s = bits(p, n, words); let t = slice::from_raw_parts(t, n);
    let (mut k, mut kept) = (0, 0u32);
    for i in i0..i1 { for j in i + 1..n {
        let d = tn93_bits_fast_inner(s[i], t[i], s[j], t[j], thr);
        if d <= thr { kept += 1; }
        *out.add(k) = d; k += 1;
    } }
    kept
}
/// One pass, all three metrics (hamming, tn93, snps) — same output layout as pairs_all_current.
#[no_mangle]
pub unsafe extern "C" fn pairs_all_bits(p: *const W, t: *const Tot, n: usize, words: usize, out: *mut f64) {
    let s = bits(p, n, words); let t = slice::from_raw_parts(t, n);
    let np = n * (n - 1) / 2; let mut k = 0;
    for i in 0..n { for j in i + 1..n {
        let d = tn93_bits_fast_inner(s[i], t[i], s[j], t[j], f64::INFINITY);
        let (sites, mism) = hamming_bits(s[i], s[j]); // cheap; a real kernel would fold this in
        *out.add(k) = if sites == 0 { 1.0 } else { mism as f64 / sites as f64 };
        *out.add(k + np) = d;
        *out.add(k + 2 * np) = mism as f64;
        k += 1;
    } }
}
#[no_mangle]
pub unsafe extern "C" fn pairs_hamming_bits(p: *const W, n: usize, words: usize, out: *mut f64) {
    let s = bits(p, n, words); let mut k = 0;
    for i in 0..n { for j in i + 1..n {
        let (sites, mism) = hamming_bits(s[i], s[j]);
        *out.add(k) = if sites == 0 { 1.0 } else { mism as f64 / sites as f64 }; k += 1;
    } }
}
