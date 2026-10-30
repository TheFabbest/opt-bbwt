use crate::{lyndon::duval, saca};

// Conceptual omega-order comparison for two byte slices
fn omega_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    // Compare their infinite expansions a.iter().cycle() and b.iter().cycle()
    let len_lcm = a.len() * b.len(); // or compare up to a.len() * b.len() steps
    let inf_a = a.iter().cycle().take(len_lcm);
    let inf_b = b.iter().cycle().take(len_lcm);
    inf_a.cmp(inf_b)
}

// TODO naive
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
            let rotation: Vec<u8> = factor.iter().cycle().skip(i).take(factor_length).cloned().collect();
            conjugates.push(rotation);
        }
    }
    
    // 3. Sort conjugates lexicographically
    conjugates.sort_by(|a, b| omega_cmp(a, b));
    
    // 4. Extract the last character of each sorted conjugate
    let mut bbwt = Vec::new();
    for conjugate in conjugates {
        if let Some(&last_char) = conjugate.last() {
            bbwt.push(last_char);
        }
    }
    bbwt
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
            let rotation: Vec<u8> = factor.iter().cycle().skip(i).take(factor_length).cloned().collect();
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


pub fn bwt(s: &[u8]) -> Vec<u8> {
    let len = s.len();
    let mut sa = vec![0; len];
    saca::optimal_suffix_array_bytes(&s, &mut sa);

    let mut bbwt = Vec::with_capacity(len); // todo can re-use sa
    for &i in &sa {
        let prev_index = if i == 0 { s.len() - 1 } else { i - 1 };
        bbwt.push(s[prev_index]);
    }
    bbwt
}