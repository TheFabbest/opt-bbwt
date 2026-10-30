use crate::saca::optimal_suffix_array_bytes;

pub fn bwt(s: &[u8]) -> Vec<u8> {
    let len = s.len();
    let mut sa = vec![0; len];
    optimal_suffix_array_bytes(&s, &mut sa);

    let mut bbwt = Vec::with_capacity(len); // todo can re-use sa
    for &i in &sa {
        let prev_index = if i == 0 { s.len() - 1 } else { i - 1 };
        bbwt.push(s[prev_index]);
    }
    bbwt
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bwt() {
        let s = b"banana";
        let expected_bwt = b"nnbaaa";
        let input = s.to_vec();
        let result = bwt(&input);
        assert_eq!(result, expected_bwt);
    }

    #[test]
    fn test_bwt_empty() {
        let s = b"";
        let expected_bwt = b"";
        let input = s.to_vec();
        let result = bwt(&input);
        assert_eq!(result, expected_bwt);
    }

    #[test]
    fn test_bwt_single_char() {
        let s = b"a";
        let expected_bwt = b"a";
        let input = s.to_vec();
        let result = bwt(&input);
        assert_eq!(result, expected_bwt);
    }
}