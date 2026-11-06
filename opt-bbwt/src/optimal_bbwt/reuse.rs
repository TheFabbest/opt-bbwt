use crate::suffix_array::circular_sais::optimal_csa;

/// Selects the conjugate of the input string that yields the least number of runs in the BBWT by simply iterating over all conjugates and computing the BBWT for each one.
/// The only optimization is that part of the algorithm for computing the circular suffix array is run only once, and the result is reused for all conjugates.
pub fn optimal_bbwt_reuse(text: &[u8]) -> usize {
    if text.is_empty() {
        return 0;
    }
    optimal_csa(text)
}
