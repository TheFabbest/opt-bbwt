use crate::{lyndon::duval, naive::omegaorder::omega_cmp};

pub fn naive_bbwt(text: &[u8]) -> Vec<u8> {
    // 1. Find Lyndon factors
    let factors = duval(&text);

    // 2. Generate all cyclic rotations (conjugates) for each factor
    let mut conjugates = Vec::new();
    for &factor_start in &factors {
        let factor_length = if let Some(&next_start) = factors.iter().find(|&&s| s > factor_start) {
            next_start - factor_start
        } else {
            text.len() - factor_start
        };
        let factor = &text[factor_start..factor_start + factor_length];
        for i in 0..factor_length {
            let rotation: Vec<u8> = factor.iter().cycle().skip(i).take(factor_length).copied().collect();
            conjugates.push(rotation);
        }
    }
    
    // 3. Sort conjugates lexicographically
    conjugates.sort_by(|a, b| omega_cmp(a, b));
    
    // 4. Extract the last character of each sorted conjugate
    let mut result = Vec::new();
    for conjugate in conjugates {
        if let Some(&last_char) = conjugate.last() {
            result.push(last_char);
        }
    }
    result
}