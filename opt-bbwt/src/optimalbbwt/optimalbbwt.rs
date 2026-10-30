use crate::bbwt::{bbwt, bbwt_with_rank}; use crate::lyndon::conjugatesfactors::{for_each_factor_of_each_conjugate, get_factor_with_conjugates_counts};
// Must return (bbwt_array, bbwt_rank)
use crate::lyndon::{duval, lyndon_suffix_table, lyndon_table};
use crate::rank::sa_to_rank;
use crate::saca::circular_suffix_array;

fn consistency_check (result: &String, sorted_factors: &Vec<String>) {
    let new_factors = duval(result.as_bytes());

    let mut count = 0;
    for &factor_start in &new_factors {
        let factor_length = if let Some(&next_start) = new_factors.iter().find(|&&s| s > factor_start) {
            next_start - factor_start
        } else {
            result.len() - factor_start
        };
        let factor = &result.as_bytes()[factor_start..factor_start + factor_length];
        
        assert!(factor == sorted_factors[count].as_bytes());
        count += 1;
    }
}

pub fn sort_all_factors(s: &[u8]) {
    // put all factors into a set, together with an array that maps each conjugate index to the number of appearances of the factor in that conjugate
    let factors = get_factor_with_conjugates_counts(s);

    // sort them in reverse lexicographical order
    let mut sorted_factors: Vec<String> = factors.keys().map(|k| String::from_utf8_lossy(k).to_string()).collect();
    sorted_factors.sort_by(|a, b| b.cmp(a));

    // concatenate them into a single string
    let result = sorted_factors.concat();
    println!("Resulting string: {}", result);
    println!("Resulting string length: {}. This is {}% of the original.", result.len(), (result.len() as f64 / s.len() as f64) * 100.0);

    // check
    consistency_check(&result, &sorted_factors);

    // compute the circular suffix array of the result
    let mut circular_sa = &mut vec![0; result.len()];
    circular_suffix_array(result.as_bytes(), &mut circular_sa);
    
    // count the number of runs in all conjugates' BBWTs
    let mut total_runs = vec![0; s.len()];
    let mut last_characters = vec![0; s.len()];
    // keep track of the built bbwt for debugging
    let mut built_bbwt = vec![vec![]; s.len()];

    // for debug, print the factors of the 5th conjugate
    // rotate s
    let rotated_s: Vec<u8> = s[4..].iter().chain(&s[..4]).cloned().collect();
    let rotated_factors = duval(&rotated_s);
    println!("Rotated string: {}", String::from_utf8_lossy(&rotated_s));
    println!("Factors of the 5th conjugate:");
    for &factor_start in &rotated_factors {
        let factor_length = if let Some(&next_start) = rotated_factors.iter().find(|&&s| s > factor_start) {
            next_start - factor_start
        } else {
            rotated_s.len() - factor_start
        };
        let factor = &rotated_s[factor_start..factor_start + factor_length];
        println!("factor {}", String::from_utf8_lossy(factor));
    }

    for i in 0..result.len() {
        // for building the bbwt, we need refer to circular_sa[i] as the starting index of the conjugate
        // the BBWT requires the predecessor of circular_sa[i] in result, rotated to the start of the lyndon factor that contains circular_sa[i]
        let start_index = circular_sa[i];
        
        // find the lyndon factor that contains start_index in result, can be obtained in sorted_factors by finding the first factor ending after start_index
        let mut factor_start = 0;
        let mut factor_length = 0;
        for factor in &sorted_factors {
            factor_length = factor.len();
            if factor_start + factor_length > start_index {
                break;
            }
            factor_start += factor_length;
        }

        let previous_index = factor_start + (start_index - factor_start + factor_length - 1) % factor_length;
        let predecessor = result.as_bytes()[previous_index];

        // now we know the predecessor and can tell whether it created a new run or not in each conjugate where the factor appears
        let conjugates = factors.get(result.as_bytes()[factor_start..factor_start + factor_length].as_ref()).unwrap();

        // debug
        println!("Computing {}-th entry in result, which refers to the conjugate {} of the factor {}. The predecessor of the conjugate is {}. The factor appears in conjugates {:?}", i, start_index, String::from_utf8_lossy(&result.as_bytes()[factor_start..factor_start + factor_length]), predecessor, conjugates);

        for conjugate_index in 0..conjugates.len() {
            let occurrences = conjugates[conjugate_index];
            if conjugates[conjugate_index] > 0 {
                // we add the predecessor to the bbwt conjugates[conjugate_index] times
                built_bbwt[conjugate_index].extend(std::iter::repeat(predecessor).take(conjugates[conjugate_index]));

                println!("Conjugate {}: last character was {}, predecessor is {} and will appear {} times in the bbwt", conjugate_index, last_characters[conjugate_index], predecessor, conjugates[conjugate_index]);
                if last_characters[conjugate_index] != predecessor {
                    total_runs[conjugate_index] += 1; // every conjugates[conjugate_index] repetitions of this factor's conjugate repeats the same predecessor, so we only count one new run, and only if the predecessor is different
                }
                last_characters[conjugate_index] = predecessor;
            }
        }
    }

    // for each conjugate, print the number of runs in its BBWT
    for i in 0..s.len() {
        println!("--------------------------------");
        println!("Conjugate starting at index {}: {} runs in BBWT", i, total_runs[i]);
        println!("In fact, its BBWT is: {:?}", bbwt(&s[i..].iter().chain(&s[..i]).cloned().collect::<Vec<u8>>()));
        println!("Built BBWT: {:?}", built_bbwt[i]);
    }
}

pub fn online_bbwt_check(input_s: &[u8]) {
    let n = input_s.len();
    if n == 0 { return; }

    // 1. FORCE THE LYNDON CONJUGATE ASSUMPTION
    // Find the least rotation to ensure the base string IS a Lyndon word.
    let mut circular_sa = vec![0; n];
    crate::saca::circular_suffix_array(input_s, &mut circular_sa);
    let least_rotation_idx = circular_sa[0];
    
    let s: Vec<u8> = input_s[least_rotation_idx..].iter()
        .chain(&input_s[..least_rotation_idx])
        .cloned()
        .collect();
    
    println!("Base Lyndon Conjugate: {}", String::from_utf8_lossy(&s));

    // 2. Compute tables on the valid Lyndon Conjugate
    crate::saca::circular_suffix_array(&s, &mut circular_sa);
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

        let current_string: Vec<u8> = s[start..].iter().chain(&s[..start]).cloned().collect();
        let rotated_string: Vec<u8> = s[new_start..].iter().chain(&s[..new_start]).cloned().collect();
        
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