//! Algorithms for selecting BBWT conjugates.

pub mod optimalbbwt;
pub mod reuse;
pub mod globalsort;
pub mod naive;

/// Selects a conjugate using the default optimized BBWT algorithm.
pub fn optimal_bbwt(text: &[u8]) -> usize {
    self::optimalbbwt::optimal_bbwt_chain(text, true)
}
