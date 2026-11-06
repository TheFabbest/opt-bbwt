pub fn lyndon_table(s: &[u8], rank: &[usize]) -> Vec<usize> {
    let len = s.len();
    let mut lyndon_table = vec![0; len];
    
    for i in (0..len).rev() {
        lyndon_table[i] = 1;
        let mut j = i + 1;
        while j < len && rank[i] < rank[j] {
            lyndon_table[i] += lyndon_table[j];
            j += lyndon_table[j];
        }
    }

    lyndon_table
}

/// Computes the Lyndon table from rank information.
pub fn compute_lyn_optimized(rank: &[usize], lyn_s: &[usize]) -> Vec<usize> {
    let n = rank.len();
    let mut lyn = vec![1; n];
    for i in (0..n).rev() {
        lyn[i] = 1;
        let mut j = i + 1;
        
        while j < n {
            let can_merge = if lyn_s[j - 1] < (j - i) {
                false 
            } else {
                rank[i] < rank[j]
            };

            if can_merge {
                lyn[i] += lyn[j];
                j += lyn[j];
            } else {
                break;
            }
        }
    }
    
    lyn
}

pub fn lyndon_suffix_table(s: &[u8]) -> Vec<usize> {
    if s.is_empty() {
        return vec![];
    }
    let mut lyndon_suffix_table = vec![0; s.len()];
    lyndon_suffix_table[0] = 1;
    let mut per = 1;
    let mut i = 0;
    for j in 1..s.len() {
        if s[j] != s[i] {
            lyndon_suffix_table[j] = j + 1;
            per = j + 1;
            i = 0;
        } else {
            lyndon_suffix_table[j] = lyndon_suffix_table[i];
            i = (i + 1) % per;
        }
    }
    lyndon_suffix_table
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lyndon_table_empty() {
        let s = b"";
        let rank = vec![];
        let expected: Vec<usize> = vec![];
        assert_eq!(lyndon_table(s, &rank), expected);
    }

    #[test]
    fn test_lyndon_table_single_char() {
        let s = b"a";
        let rank = vec![0];
        let expected = vec![1];
        assert_eq!(lyndon_table(s, &rank), expected);
    }

    #[test]
    fn test_lyndon_table_simple() {
        let s = b"banana";
        let rank = vec![3, 2, 5, 1, 4, 0];
        let expected = vec![1, 2, 1, 2, 1, 1];
        assert_eq!(lyndon_table(s, &rank), expected);
    }

    #[test]
    fn test_lyndon_suffix_table_empty() {
        let s = b"";
        let expected: Vec<usize> = vec![];
        assert_eq!(lyndon_suffix_table(s), expected);
    }

    #[test]
    fn test_lyndon_suffix_table_single_char() {
        let s = b"a";
        let expected = vec![1];
        assert_eq!(lyndon_suffix_table(s), expected);
    }

    #[test]
    fn test_lyndon_suffix_table_simple() {
        let s = b"banana";
        let expected = vec![1, 2, 3, 4, 5, 6];
        assert_eq!(lyndon_suffix_table(s), expected);
    }

    #[test]
    fn test_lyndon_table_optimized_empty() {
        let rank = vec![];
        let lyn_s = vec![];
        let expected: Vec<usize> = vec![];
        assert_eq!(compute_lyn_optimized(&rank, &lyn_s), expected);
    }

    #[test]
    fn test_lyndon_table_optimized_single_char() {
        let rank = vec![0];
        let lyn_s = vec![1];
        let expected = vec![1];
        assert_eq!(compute_lyn_optimized(&rank, &lyn_s), expected);
    }

    #[test]
    fn test_lyndon_table_optimized_simple() {
        let rank = vec![3, 2, 5, 1, 4, 0];
        let lyn_s = vec![1, 2, 1, 2, 1, 1];
        let expected = vec![1, 2, 1, 2, 1, 1];
        assert_eq!(compute_lyn_optimized(&rank, &lyn_s), expected);
    }
}