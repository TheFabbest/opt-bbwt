pub mod optimalbbwt;
pub mod reuse;
pub mod globalsort;
pub mod naive;

pub fn optimal_bbwt(text: &[u8]) -> usize {
    self::optimalbbwt::optimal_bbwt_chain(text, true)
}

pub mod variants {
    pub use super::optimalbbwt::{online_bbwt_check, optimal_bbwt_chain};
    pub use super::reuse::optimal_bbwt_reuse;
    pub use super::globalsort::optimal_bbwt_global_sort;
    pub use super::naive::naive_optimal_bbwt;
}