use crate::suffix_array::suffix_array;

/// Computes the Burrows-Wheeler Transform (BWT) of a string.
pub fn bwt(s: &[u8]) -> Vec<u8> {
    let mut output = vec![0; s.len()];
    bwt_into(s, &mut output);
    output
}

/// Computes the BWT into a caller-provided buffer.
///
/// The buffer must have the same length as `s`.
pub fn bwt_into(s: &[u8], output: &mut [u8]) {
    assert_eq!(
        output.len(),
        s.len(),
        "BWT output buffer must match input length"
    );
    let sa = suffix_array(s);
    for (output_index, &i) in sa.iter().enumerate() {
        let prev_index = if i == 0 { s.len() - 1 } else { i - 1 };
        output[output_index] = s[prev_index];
    }
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

    #[test]
    fn test_bwt_into() {
        let output = bwt(b"banana");
        assert_eq!(&output, b"nnbaaa");
    }
}
