// this implementation refers to https://github.com/hmmpiatkowski/bbwt

use crate::{lyndon::duval, saca::circular_suffix_array};

/// Computes the Bijective Burrows-Wheeler Transform (BBWT) of a string.
pub fn bbwt(text: &[u8]) -> Vec<u8> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }

    let factors = duval(text);
    let csa = circular_suffix_array(text);

    // Map global indices back to their Lyndon factor boundaries to pick
    // the preceding character in cyclic order.
    let transformed: Vec<u8> = csa
        .into_iter()
        .map(|idx| {
            // Find which Lyndon factor this index belongs to (TODO store them)
            let factor_start = factors
                .iter()
                .rev()
                .find(|&&start| start <= idx)
                .copied()
                .unwrap_or(0);
            let factor_len = if let Some(&next_start) = factors.iter().find(|&&s| s > factor_start) {
                next_start - factor_start
            } else {
                n - factor_start
            };

            // Preceding character within its Lyndon factor context
            let prev_idx = if idx == factor_start {
                factor_start + factor_len - 1
            } else {
                idx - 1
            };

            text[prev_idx]
        })
        .collect();

    transformed
}

/// Computes the inverse of the Bijective Burrows-Wheeler Transform (BBWT).
pub fn unbbwt(in_str: &[u8]) -> Vec<u8> {
    let len = in_str.len();
    if len == 0 {
        return Vec::new();
    }
    if len == 1 {
        return vec![in_str[0]];
    }

    let alph_begin = in_str.iter().copied().min().unwrap_or(0);
    let alph_size = in_str.iter().copied().max().unwrap_or(0) as usize + 1;
    let mut chars_count = vec![0; alph_size];
    let mut chars_before = vec![0; alph_size];
    let mut chars_seen = vec![0; alph_size];
    let mut std_perm = vec![0; len];

    // lambda to map char to index in the alphabet
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
    let mut out_str = vec![0; len];

    for j in 0..len {
        if std_perm[j] == std::usize::MAX {
            continue;
        }

        let mut in_pos = j;

        while std_perm[in_pos] != std::usize::MAX {
            out_str[out_pos] = in_str[in_pos];
            if out_pos == 0 { break; }
            out_pos -= 1;
            let t = in_pos;
            in_pos = std_perm[in_pos];
            std_perm[t] = std::usize::MAX;
        }
    }

    out_str
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
    fn test_bbwt_vs_naive() {
        let test_strings = vec![
            b"banana".to_vec(),
            b"abracadabra".to_vec(),
            b"".to_vec(),
            b"a".to_vec()];
        for s in test_strings {
            let transformed = bbwt(&s);
            let naive_transformed = crate::naive::bbwt::bbwt(&s);
            assert_eq!(transformed, naive_transformed, "BBWT mismatch for input: {:?}", String::from_utf8_lossy(&s));
        }
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