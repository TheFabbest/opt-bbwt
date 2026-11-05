use crate::{lyndon::duval};

fn omega_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    // Compare their infinite expansions a.iter().cycle() and b.iter().cycle()
    let len_lcm = a.len() * b.len(); // or compare up to a.len() * b.len() steps
    let inf_a = a.iter().cycle().take(len_lcm);
    let inf_b = b.iter().cycle().take(len_lcm);
    inf_a.cmp(inf_b)
}

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