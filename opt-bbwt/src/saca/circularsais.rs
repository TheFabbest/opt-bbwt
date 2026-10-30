use crate::{naive::circular_suffix_array_naive};

/// Builds the Circular Suffix Array (CSA) by creating doubled Lyndon factors
/// and ordering all cyclic rotations in true BBWT order.
pub fn circular_suffix_array(bytes: &[u8]) -> Vec<usize> {
    return circular_suffix_array_naive(bytes);
}


#[cfg(test)]
mod tests {
    use crate::naive::csa::circular_suffix_array_naive;

use super::*;

    #[test]
    fn test_circular_suffix_array_empty() {
        let s = b"";
        let csa = circular_suffix_array(s);
        assert_eq!(csa, vec![]);
    }

    #[test]
    fn test_circular_suffix_array_single() {
        let s = b"z";
        let csa = circular_suffix_array(s);
        assert_eq!(csa, vec![0]);
    }

    #[test]
    fn test_circular_suffix_array() {
        let s = b"banana";
        let csa = circular_suffix_array(s);
        assert_eq!(csa, vec![5, 3, 1, 0, 4, 2]);
    }

    #[test]
    fn test_circular_suffix_array_compare_with_naive() {
        let s = b"banana";
        let csa = circular_suffix_array(s);
        let naive_csa = circular_suffix_array_naive(s);
        assert_eq!(csa, naive_csa);
    }
}