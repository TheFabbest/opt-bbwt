use crate::bbwt::{bbwt}; use crate::lyndon::conjugatesfactors::{get_factor_with_conjugates_grouped};
use crate::lyndon::{lyndon_suffix_table, lyndon_table, shiloach};
use crate::naive::bbwt::bbwt_with_rank;
use crate::rank::{sa_to_rank};
use crate::saca::circularsais::optimal_csa;
use crate::saca::{circular_suffix_array};

pub fn naive_optimal_bbwt(text: &[u8]) -> usize {
    let (start, n) = shiloach(text);

    let s = (0..n).map(|k| text[(start + k) % text.len()]).collect::<Vec<u8>>();

    let mut min_runs_index = 0;
    let mut min_runs = usize::MAX;
    for i in 0..n {
        let conjugate: Vec<u8> = s[i..].iter().chain(&s[..i]).copied().collect();
        let bbwt_result = bbwt(&conjugate);
        let runs = bbwt_result.windows(2).filter(|w| w[0] != w[1]).count() + 1; // count runs in BBWT
        if runs < min_runs {
            min_runs = runs;
            min_runs_index = i;
        }
    }

    min_runs_index
}

pub fn optimal_bbwt_through_w(text: &[u8]) -> usize {
    let (start, period_length) = shiloach(text);
    let s = &(0..period_length).map(|k| text[(start + k) % text.len()]).collect::<Vec<u8>>();

    // Directly get the flattened, sorted arena
    let grouped = get_factor_with_conjugates_grouped(s);
    let sorted = grouped.sorted_factors;

    // W = concatenation
    let mut w: Vec<u8> = Vec::new();
    let mut starts = Vec::with_capacity(sorted.len());
    for &f in &sorted {  
        starts.push(w.len());
        w.extend_from_slice(f);
    }

    let mut factor_of = vec![0usize; w.len()];
    for (id, f) in sorted.iter().enumerate() {
        for p in starts[id]..starts[id] + f.len() {
            factor_of[p] = id;
        }
    }

    // Compute the circular suffix array of W
    let circular_sa = circular_suffix_array(&w);
    
    // Use 256 as a sentinel for "None" to avoid Option<u8> memory overhead and branching
    let mut run_state: Vec<(usize, u16)> = vec![(0, 256u16); s.len()];
    
    for &p in &circular_sa {
        let id = factor_of[p];
        let f = sorted[id];
        
        // REMOVED MODULO: Calculate the predecessor index safely
        let offset = p - starts[id];
        let pred = if offset == 0 {
            f[f.len() - 1]
        } else {
            f[offset - 1]
        } as u16;
        
        let idx_start = grouped.factor_starts[id];
        let idx_len = grouped.factor_lengths[id];
        let indices = &grouped.factor_indices[idx_start .. idx_start + idx_len];
        
        for &c in indices {
            let state = &mut run_state[c];
            if state.1 != pred {
                state.0 += 1;
                state.1 = pred;
            }
        }
    }

    run_state.iter().enumerate().min_by_key(|&(_, state)| state.0).map(|(i, _)| i).unwrap()
}

pub fn optimal_bbwt_reuse(text: &[u8]) -> usize {
    optimal_csa(text)
}

// WIP
pub fn online_bbwt_check(input_s: &[u8]) {
    let n = input_s.len();
    if n == 0 { return; }

    // 1. FORCE THE LYNDON CONJUGATE ASSUMPTION
    // Find the least rotation to ensure the base string IS a Lyndon word.
    let mut circular_sa = circular_suffix_array(input_s);
    let least_rotation_idx = circular_sa[0];
    
    let s: Vec<u8> = input_s[least_rotation_idx..].iter()
        .chain(&input_s[..least_rotation_idx])
        .copied()
        .collect();
    
    println!("Base Lyndon Conjugate: {}", String::from_utf8_lossy(&s));

    // 2. Compute tables on the valid Lyndon Conjugate
    circular_sa = circular_suffix_array(&s);
    let lyn_suffix_table = lyndon_suffix_table(&s);
    let sa_rank = sa_to_rank(&circular_sa);
    let lyn_table = lyndon_table(&s, &sa_rank); 

    // Bucket cumulative counts for O(1) BBWT insertion
    let mut c = vec![0; 256];
    for &byte in &s { c[byte as usize] += 1; }
    let mut sum = 0;
    for count in c.iter_mut() {
        let current = *count;
        *count = sum;
        sum += current;
    }

    for shift in 0..n {
        let start = (n - shift) % n;
        let new_start = (n - shift - 1) % n;
        let removed_char = s[new_start]; 

        let current_string: Vec<u8> = s[start..].iter().chain(&s[..start]).copied().collect();
        let rotated_string: Vec<u8> = s[new_start..].iter().chain(&s[..new_start]).copied().collect();
        
        let (mut current_bbwt, bbwt_rank) = bbwt_with_rank(&current_string);
        let (new_bbwt, _) = bbwt_with_rank(&rotated_string);

        let l = if start == 0 {
            let mut curr = 0;
            let mut last_len = 0;
            while curr < n {
                last_len = lyn_table[curr];
                curr += last_len;
            }
            last_len
        } else {
            lyn_suffix_table[start - 1]
        };

        let prepend_is_safe = s[new_start] >= s[start];
        let remove_is_safe = if l == 1 {
            true 
        } else {
            let char_before_removed = (start + n - 2) % n;
            lyn_suffix_table[char_before_removed] == l - 1
        };

        if prepend_is_safe && remove_is_safe {
            // 1. Remove the entry corresponding to the cyclic shift starting at the removed char
            let r_last = bbwt_rank[n - 1];
            current_bbwt.remove(r_last);

            // 2. Insert the removed character at the head of its bucket
            let bucket_idx = c[removed_char as usize];
            current_bbwt.insert(bucket_idx, removed_char);

            // 3. Update predecessor ONLY if the factor wasn't entirely destroyed
            if l > 1 {
                let mut r_first = bbwt_rank[n - l];

                // Compensate for array shifts
                if r_last < r_first { r_first -= 1; }
                if bucket_idx <= r_first { r_first += 1; }

                current_bbwt[r_first] = current_string[n - 2];
            }

            if current_bbwt != new_bbwt {
                println!("Shift {}: [!] BUG: Expected {:?}, got {:?}", shift, new_bbwt, current_bbwt);
            } else {
                println!("Shift {}: [+] Patch successful!", shift);
            }
        } else {
            println!("Shift {}: Cascade detected. Recomputing.", shift);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optimal_bbwt_through_w_vs_naive() {
        let test_strings: Vec<&[u8]> = vec![
            b"abcabcabc",
            b"aaaaaa",
            b"abababab",
            b"xyzxyzxyz",
            b"abcdabcdabcd",
            b"abcdeabcdeabcde",
            b"aabbccddeeff",
            b"thequickbrownfoxjumpsoverthelazydog",
            b"banana",
            b"mississippi",
        ];

        for s in test_strings {
            let result = optimal_bbwt_through_w(s);
            let naive_result = naive_optimal_bbwt(s);

            // get result-th rotation of s
            let rotation: Vec<u8> = s[result..].iter().chain(&s[..result]).copied().collect();
            let naive_rotation: Vec<u8> = s[naive_result..].iter().chain(&s[..naive_result]).copied().collect();

            // calculate bbwt for both rotations
            let bbwt_result = bbwt(&rotation);
            let bbwt_naive_result = bbwt(&naive_rotation);

            // count runs for both bbwt results
            let runs_result = bbwt_result.windows(2).filter(|w| w[0] != w[1]).count() + 1;
            let runs_naive_result = bbwt_naive_result.windows(2).filter(|w| w[0] != w[1]).count() + 1;

            assert_eq!(runs_result, runs_naive_result, "Mismatch in runs for string: {:?}. Naive found optimal rotation {} with BBWT {}. Through_W found {} with BBWT {}", String::from_utf8_lossy(s), String::from_utf8_lossy(&naive_rotation), String::from_utf8_lossy(&bbwt_naive_result), String::from_utf8_lossy(&rotation), String::from_utf8_lossy(&bbwt_result));
        }
    }
}