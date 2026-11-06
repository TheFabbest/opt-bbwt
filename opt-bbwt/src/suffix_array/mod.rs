pub mod circular_sais;
pub mod naive_circular;
pub mod optimal;
pub mod sais;

pub fn suffix_array(text: &[u8]) -> Vec<usize> {
    optimal::optimal_suffix_array_bytes(text)
}

pub fn circular_suffix_array(text: &[u8]) -> Vec<usize> {
    circular_sais::circular_suffix_array(text)
}
