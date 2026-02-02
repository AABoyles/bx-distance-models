/// Utility functions for nucleotide analysis

#[inline]
pub fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

#[inline]
pub fn is_purine(b: u8) -> bool {
    matches!(b, b'A' | b'G')
}

#[inline]
pub fn is_pyrimidine(b: u8) -> bool {
    matches!(b, b'C' | b'T')
}
