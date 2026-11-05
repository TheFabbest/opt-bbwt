use crate::{lyndon::duval};

fn omega_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    // Compare their infinite expansions a.iter().cycle() and b.iter().cycle()
    let len_lcm = a.len() * b.len(); // or compare up to a.len() * b.len() steps
    let inf_a = a.iter().cycle().take(len_lcm);
    let inf_b = b.iter().cycle().take(len_lcm);
    inf_a.cmp(inf_b)
}

// circular suffix array
pub fn naive_circular_suffix_array(s: &[u8]) -> Vec<usize> {
    
    let mut csa = vec![0; s.len()];

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
        csa[i] = start;
    }

    csa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular_suffix_array_naive_empty() {
        let s = b"";
        let expected: Vec<usize> = vec![];
        assert_eq!(naive_circular_suffix_array(s), expected);
    }

    #[test]
    fn test_circular_suffix_array_naive_single_char() {
        let s = b"a";
        let expected = vec![0];
        assert_eq!(naive_circular_suffix_array(s), expected);
    }

    #[test]
    fn test_circular_suffix_array_naive() {
        let s = b"banana";
        let expected = vec![5, 3, 1, 0, 4, 2];
        assert_eq!(naive_circular_suffix_array(s), expected);
    }
}