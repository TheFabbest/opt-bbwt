/// Computes the starting indices of the Lyndon factors of `s` using Duval's algorithm in O(n) time and O(1) extra space.
pub fn duval(s: &[u8]) -> Vec<usize> {
    let n = s.len();
    let mut result = Vec::new();
    let mut i = 0;

    while i < n {
        let mut j = i + 1;
        let mut k = i;

        while j < n && s[k] <= s[j] {
            if s[k] < s[j] {
                k = i;
            } else {
                k += 1;
            }
            j += 1;
        }

        while i <= k {
            result.push(i);
            i += j - k;
        }
    }

    result
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_duval_empty() {
        let s = b"";
        let expected: Vec<usize> = vec![];
        assert_eq!(duval(s), expected);
    }

    #[test]
    fn test_duval_single_char() {
        let s = b"a";
        let expected = vec![0];
        assert_eq!(duval(s), expected);
    }

    #[test]
    fn test_duval_repeated_chars() {
        let s = b"aaa";
        let expected = vec![0,1,2];
        assert_eq!(duval(s), expected);
    }

    #[test]
    fn test_duval_distinct_chars() {
        let s = b"banana";
        let expected = vec![0, 1, 3, 5];
        assert_eq!(duval(s), expected);
    }
}