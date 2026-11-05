pub mod optimalsaca;
pub mod circularsais;
pub mod sais;
pub mod naivecsa;

pub fn suffix_array(text: &[u8]) -> Vec<usize> {
    optimalsaca::optimal_suffix_array_bytes(text)
}

pub fn circular_suffix_array(text: &[u8]) -> Vec<usize> {
    circularsais::circular_suffix_array(text)
}

pub mod variants {
    pub use super::optimalsaca::{optimal_suffix_array_bytes, optimal_suffix_array_bytes_into};
    pub use super::circularsais::{circular_suffix_array, circular_suffix_array_into};
    pub use super::sais::{sais, sais_into};
    pub use super::naivecsa::{naive_circular_suffix_array};
}