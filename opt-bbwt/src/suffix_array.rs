//! Suffix-array construction algorithms.
//!
//! The module exposes the default suffix-array implementations used by the
//! library, as well as alternative implementations under their respective
//! submodules. Suffix arrays contain the starting positions of all suffixes
//! in lexicographic order.

pub mod circular_sais;
pub mod naive_circular;
pub mod optimal;
pub mod sais;

/// Builds a suffix array for `text`.
///
/// The returned positions refer to suffixes of `text` and are ordered
/// lexicographically. The default implementation is the optimized
/// in-place algorithm.
///
/// # Examples
///
/// ```
/// use opt_bbwt::suffix_array::suffix_array;
///
/// assert_eq!(suffix_array(b"banana"), vec![5, 3, 1, 0, 4, 2]);
/// ```
pub fn suffix_array(text: &[u8]) -> Vec<usize> {
    optimal::optimal_suffix_array_bytes(text)
}

/// Builds a suffix array for all cyclic rotations of `text`.
///
/// The returned positions identify the starting position of each rotation in
/// lexicographic order. An empty input produces an empty suffix array.
///
/// # Examples
///
/// ```
/// use opt_bbwt::suffix_array::circular_suffix_array;
///
/// assert_eq!(circular_suffix_array(b"banana"), vec![5, 3, 1, 0, 4, 2]);
/// ```
pub fn circular_suffix_array(text: &[u8]) -> Vec<usize> {
    circular_sais::circular_suffix_array(text)
}

/// Converts suffix-array positions into ranks.
///
/// # Examples
/// ```
/// use opt_bbwt::suffix_array::{suffix_array, sa_to_rank};
/// let text = b"banana";
/// let sa = suffix_array(text);
/// let rank = sa_to_rank(&sa);
/// assert_eq!(rank, vec![3, 2, 5, 1, 4, 0]);
/// ```
pub fn sa_to_rank(sa: &[usize]) -> Vec<usize> {
    let n = sa.len();
    let mut rank = vec![0; n];
    for (i, &suffix_index) in sa.iter().enumerate() {
        rank[suffix_index] = i;
    }
    rank
}
