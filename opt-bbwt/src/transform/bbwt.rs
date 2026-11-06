use crate::{lyndon::duval, suffix_array::circular_suffix_array};

pub mod naive;

/// Computes the Bijective Burrows-Wheeler Transform (BBWT) of `text`.
///
/// The returned vector has the same length as the input.
///
/// # Examples
///
/// ```
/// use opt_bbwt::transform::bbwt::bbwt;
///
/// assert_eq!(bbwt(b""), b"");
/// ```
pub fn bbwt(text: &[u8]) -> Vec<u8> {
    let mut output = vec![0; text.len()];
    bbwt_into(text, &mut output);
    output
}

/// Computes the BBWT into a caller-provided buffer.
///
/// `output` must have the same length as `text`.
///
/// # Panics
///
/// Panics if `output.len() != text.len()`.
pub fn bbwt_into(text: &[u8], output: &mut [u8]) {
    assert_eq!(
        output.len(),
        text.len(),
        "BBWT output buffer must match input length"
    );
    let n = text.len();
    if n == 0 {
        return;
    }
    let factors = duval(text);
    let csa = circular_suffix_array(text);

    for (output_index, idx) in csa.into_iter().enumerate() {
        let factor_start = factors
            .iter()
            .rev()
            .find(|&&start| start <= idx)
            .copied()
            .unwrap_or(0);
        let factor_len = factors
            .iter()
            .find(|&&start| start > factor_start)
            .map_or(n - factor_start, |&next_start| next_start - factor_start);
        let prev_idx = if idx == factor_start {
            factor_start + factor_len - 1
        } else {
            idx - 1
        };
        output[output_index] = text[prev_idx];
    }
}

/// Computes the inverse of the Bijective Burrows-Wheeler Transform (BBWT).
///
/// For BBWT output, this reconstructs the original input:
/// `unbbwt(&bbwt(text)) == text`.
///
/// # Examples
///
/// ```
/// use opt_bbwt::transform::bbwt::{bbwt, unbbwt};
///
/// let input = b"banana";
/// assert_eq!(unbbwt(&bbwt(input)), input);
/// ```
pub fn unbbwt(in_str: &[u8]) -> Vec<u8> {
    let mut output = vec![0; in_str.len()];
    unbbwt_into(in_str, &mut output);
    output
}

/// Computes the inverse BBWT into a caller-provided buffer.
///
/// `output` must have the same length as `in_str`.
///
/// # Panics
///
/// Panics if `output.len() != in_str.len()`.
pub fn unbbwt_into(in_str: &[u8], output: &mut [u8]) {
    assert_eq!(
        output.len(),
        in_str.len(),
        "inverse BBWT output buffer must match input length"
    );
    let len = in_str.len();
    if len == 0 {
        return;
    }
    if len == 1 {
        output[0] = in_str[0];
        return;
    }

    let alph_begin = in_str.iter().copied().min().unwrap_or(0);
    let alph_size = in_str.iter().copied().max().unwrap_or(0) as usize + 1;
    let mut chars_count = vec![0; alph_size];
    let mut chars_before = vec![0; alph_size];
    let mut chars_seen = vec![0; alph_size];
    let mut std_perm = vec![0; len];

    let char_to_index = |ch: u8| (ch - alph_begin) as usize;

    for &ch in in_str {
        chars_count[char_to_index(ch)] += 1;
    }

    for i in 1..alph_size {
        chars_before[i] = chars_before[i - 1] + chars_count[i - 1];
    }

    for (i, &ch) in in_str.iter().enumerate() {
        std_perm[i] = chars_before[char_to_index(ch)] + chars_seen[char_to_index(ch)];
        chars_seen[char_to_index(ch)] += 1;
    }

    let mut out_pos = len - 1;
    for j in 0..len {
        if std_perm[j] == std::usize::MAX {
            continue;
        }

        let mut in_pos = j;

        while std_perm[in_pos] != std::usize::MAX {
            output[out_pos] = in_str[in_pos];
            if out_pos == 0 {
                break;
            }
            out_pos -= 1;
            let t = in_pos;
            in_pos = std_perm[in_pos];
            std_perm[t] = std::usize::MAX;
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bbwt_empty() {
        let s = b"";
        let transformed = bbwt(s);
        assert_eq!(transformed, b"");
    }

    #[test]
    fn test_bbwt_single_char() {
        let s = b"a";
        let transformed = bbwt(s);
        assert_eq!(transformed, b"a");
    }

    #[test]
    fn test_unbbwt_empty() {
        let s = b"";
        let restored = unbbwt(s);
        assert_eq!(restored, b"");
    }

    #[test]
    fn test_unbbwt_single_char() {
        let s = b"a";
        let restored = unbbwt(s);
        assert_eq!(restored, b"a");
    }

    #[test]
    fn test_unbbwt_simple() {
        let s = b"banana";
        let transformed = bbwt(s);
        let restored = unbbwt(&transformed);
        assert_eq!(restored, s);
    }
}
