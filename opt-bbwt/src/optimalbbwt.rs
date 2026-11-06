pub mod optimalbbwt;
pub mod reuse;
pub mod globalsort;
pub mod naive;

/// Computes the default optimized BBWT conjugate.
pub fn optimal_bbwt(text: &[u8]) -> usize {
    self::optimalbbwt::optimal_bbwt_chain(text, true)
}
