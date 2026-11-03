use crate::{lyndon::duval, naive::omegaorder::omega_cmp};

// circular suffix array
pub fn circular_suffix_array_naive(s: &[u8]) -> Vec<usize> {
    
    let mut sa = vec![0; s.len()];

    // duval
    let factors = duval(s);

    // get all rotations of all factors, keep track of their starting indices in the original string
    let mut rotations = Vec::new();
    
    // in reverse because if the same factor appears multiple times, we want the last one to be considered first in the circular suffix array
    for &factor_start in factors.iter().rev() {
        let factor_len = if let Some(&next_start) = factors.iter().find(|&&s| s > factor_start) {
            next_start - factor_start
        } else {
            s.len() - factor_start
        };
        let factor = &s[factor_start..factor_start + factor_len];
        for i in 0..factor_len {
            let rotation = factor[i..].iter().chain(factor[..i].iter()).copied().collect::<Vec<u8>>();
            rotations.push((rotation, factor_start + i));
        }
    }

    // sort rotations according to omega order
    rotations.sort_by(|a, b| omega_cmp(&a.0, &b.0));

    // populate the suffix array with the starting indices of the sorted rotations
    for (i, (_, start)) in rotations.into_iter().enumerate() {
        sa[i] = start;
    }

    sa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular_suffix_array_naive_empty() {
        let s = b"";
        let expected: Vec<usize> = vec![];
        assert_eq!(circular_suffix_array_naive(s), expected);
    }

    #[test]
    fn test_circular_suffix_array_naive_single_char() {
        let s = b"a";
        let expected = vec![0];
        assert_eq!(circular_suffix_array_naive(s), expected);
    }

    #[test]
    fn test_circular_suffix_array_naive() {
        let s = b"banana";
        let expected = vec![5, 3, 1, 0, 4, 2];
        assert_eq!(circular_suffix_array_naive(s), expected);
    }
}