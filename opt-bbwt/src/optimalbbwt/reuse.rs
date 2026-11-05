use crate::suffix_array::circularsais::optimal_csa;

pub fn optimal_bbwt_reuse(text: &[u8]) -> usize {
    optimal_csa(text)
}