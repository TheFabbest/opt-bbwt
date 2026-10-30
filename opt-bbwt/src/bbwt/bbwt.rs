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
            // Find which Lyndon factor this index belongs to
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