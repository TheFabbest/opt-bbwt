//! Exact algorithms for selecting the conjugate of a string that has the least number of runs in its BBWT.
//! In practice, you may often want to use the heuristics in the `heuristics` module instead.

pub mod chain;
pub mod globalsort;
pub mod naive;
pub mod reuse;

/// Selects a conjugate using the default optimized BBWT algorithm.
pub fn optimal_bbwt(text: &[u8]) -> usize {
    self::chain::optimal_bbwt_chain(text, true)
}
