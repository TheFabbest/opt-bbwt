use crate::{lyndon::duval, naive::omegaorder::omega_cmp};

pub fn bbwt(text: &[u8]) -> Vec<u8> {
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


pub fn bbwt_with_rank(text: &[u8]) -> (Vec<u8>, Vec<usize>) {
    // get all conjugates of the lyndon factors of the string, sort them while keeping track of their original indices, and return the last character of each sorted conjugate as the BBWT, along with the rank array
    let factors = duval(&text);
    let mut conjugates_with_index = Vec::new();
    for &factor_start in &factors {
        let factor_length = if let Some(&next_start) = factors.iter().find(|&&s| s > factor_start) {
            next_start - factor_start
        } else {
            text.len() - factor_start
        };
        let factor = &text[factor_start..factor_start + factor_length];
        for i in 0..factor_length {
            let rotation: Vec<u8> = factor.iter().cycle().skip(i).take(factor_length).copied().collect();
            conjugates_with_index.push((rotation, factor_start + i));
        }
    }

    // sort conjugates lexicographically while keeping track of their original indices
    conjugates_with_index.sort_by(|(a, _), (b, _)| omega_cmp(a, b));

    let mut bbwt = Vec::new();
    let mut rank = Vec::new();
    for (conjugate, original_index) in conjugates_with_index {
        if let Some(&last_char) = conjugate.last() {
            bbwt.push(last_char);
        }
        rank.push(original_index);
    }

    (bbwt, rank) // Return both the BBWT and the rank array
}