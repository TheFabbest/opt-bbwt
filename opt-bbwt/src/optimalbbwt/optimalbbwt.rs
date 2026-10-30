use crate::bbwt::{bbwt, bbwt_with_rank}; use crate::lyndon::conjugatesfactors::{get_factor_with_conjugates_counts};
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
        
        assert!(factor == sorted_factors[count].as_bytes(), "Consistency check failed: factor {:?} does not match sorted factor {:?}", String::from_utf8_lossy(factor), sorted_factors[count]);
        count += 1;
    }
}

pub fn sort_all_factors(s: &[u8]) -> usize {
    let correctness_checks = false; // set to true to enable correctness checks, which will slow down the program significantly

    // put all factors into a set, together with an array that maps each conjugate index to the number of appearances of the factor in that conjugate
    let factors = get_factor_with_conjugates_counts(s);

    // sort them in reverse lexicographical order (it's the same as omega-order because they are all Lyndon words)
    let mut sorted_factors: Vec<String> = factors.keys().map(|k| String::from_utf8_lossy(k).to_string()).collect();
    sorted_factors.sort_by(|a, b| b.cmp(a));

    // concatenate them into a single string
    let result = sorted_factors.concat();

    // check
    if correctness_checks {
        consistency_check(&result, &sorted_factors);
    }

    // compute the circular suffix array of the result
    let mut circular_sa = &mut vec![0; result.len()];
    circular_suffix_array(result.as_bytes(), &mut circular_sa);
    
    // count the number of runs in all conjugates' BBWTs
    let mut total_runs = vec![0; s.len()];
    let mut last_characters = vec![0; s.len()];

    // keep track of the built bbwt for debugging
    let mut built_bbwt= vec![vec![]; s.len()];

    for i in 0..result.len() {
        // for building the bbwt, we need refer to circular_sa[i] as the starting index of the conjugate
        // the BBWT requires the predecessor of circular_sa[i] in result, rotated to the start of the lyndon factor that contains circular_sa[i]
        let start_index = circular_sa[i];
        
        // find the lyndon factor that contains start_index in result, can be obtained in sorted_factors by finding the first factor ending after start_index
        // linear on the number of factors, can be improved
        let mut factor_start = 0;
        let mut factor_length = 0;
        for factor in &sorted_factors {
            factor_length = factor.len();
            if factor_start + factor_length > start_index {
                break;
            }
            factor_start += factor_length;
        }

        // the previous index is the index of the character before the start_index in the current lyndon factor, wrapping around the factor if necessary
        let previous_index = factor_start + (start_index - factor_start + factor_length - 1) % factor_length;
        let predecessor = result.as_bytes()[previous_index];

        // now we know the predecessor and can tell whether it created a new run or not in each conjugate where the factor appears
        let conjugates = factors.get(result.as_bytes()[factor_start..factor_start + factor_length].as_ref()).unwrap();

        for conjugate_index in 0..conjugates.len() {
            let occurrences = conjugates[conjugate_index];
            if occurrences > 0 {
                // we add the predecessor to the bbwt occurrences times
                if correctness_checks {
                    built_bbwt[conjugate_index].extend(std::iter::repeat(predecessor).take(occurrences));
                }

                if last_characters[conjugate_index] != predecessor {
                    total_runs[conjugate_index] += 1; // every occurrences repetitions of this factor's conjugate repeats the same predecessor, so we only count one new run, and only if the predecessor is different
                }
                last_characters[conjugate_index] = predecessor;
            }
        }
    }


    if correctness_checks {
        // for each conjugate, print the number of runs in its BBWT
        for i in 0..s.len() {
            println!("--------------------------------");
            println!("Conjugate starting at index {}: {} runs in BBWT", i, total_runs[i]);
            let its_bbwt = built_bbwt[i].iter().collect::<Vec<&u8>>().iter().map(|&&c| c as char).collect::<String>();
            println!("In fact, its BBWT is: {:?}", its_bbwt);
            assert_eq!(its_bbwt, String::from_utf8_lossy(&bbwt(&s[i..].iter().chain(&s[..i]).cloned().collect::<Vec<u8>>())));
        }
    }

    // return the index of the max in total_runs
    total_runs.iter().enumerate().max_by_key(|&(_, count)| count).map(|(index, _)| index).unwrap()
}

pub fn naive_optimal_bbwt(s: &[u8]) -> usize {
    let n = s.len();
    
    let mut max_runs_index = 0;
    let mut max_runs = 0;
    for i in 0..n {
        let conjugate: Vec<u8> = s[i..].iter().chain(&s[..i]).cloned().collect();
        let bbwt_result = bbwt(&conjugate);
        let runs = bbwt_result.windows(2).filter(|w| w[0] != w[1]).count() + 1; // count runs in BBWT
        if runs > max_runs {
            max_runs = runs;
            max_runs_index = i;
        }
    }

    max_runs_index
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