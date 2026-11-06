use crate::suffix_array::circular_sais::optimal_csa;

pub fn optimal_bbwt_reuse(text: &[u8]) -> usize {
    if text.is_empty() {
        return 0;
    }
    optimal_csa(text)
}