pub fn from_period_slice(s: &[u8], start: usize, period_length: usize) -> Vec<u8> {
    let len = s.len();
    let word: Vec<u8> = if start + period_length <= len {
        s[start..start + period_length].to_vec()
    } else {
        [&s[start..], &s[..start + period_length - len]].concat()
    };
    word
}

pub fn rotate_string(s: &[u8], start: usize) -> Vec<u8> {
    let len = s.len();
    let start = start % len;
    let rotated: Vec<u8> = [&s[start..], &s[..start]].concat();
    rotated
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_period_slice_empty() {
        let s = vec![];
        let result = from_period_slice(&s, 0, 0);
        assert_eq!(result, vec![]);
    }

    #[test]
    fn test_from_period_slice_single_char() {
        let s = vec![1];
        let result = from_period_slice(&s, 0, 1);
        assert_eq!(result, vec![1]);
    }

    #[test]
    fn test_from_period_slice_multiple_chars() {
        let s = vec![1, 2, 3, 4, 5];
        let result = from_period_slice(&s, 4, 5);
        assert_eq!(result, vec![5,1,2,3,4]);
    }

    #[test]
    fn test_from_period_slice_periodic() {
        let s = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let result = from_period_slice(&s, 5, 4);
        assert_eq!(result, vec![2, 3, 4, 1]);
    }
}