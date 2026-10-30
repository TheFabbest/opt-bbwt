use crate::lyndon::{self, CanonicalRotationFinder, shiloach};
use crate::optimalbbwt::{online_bbwt_check};
use crate::string::rotate_string;
use crate::{rank, saca};

// takes a lambda and executes it on each factor of each conjugate of the string
pub fn for_each_factor_of_each_conjugate<F>(s: &[u8], mut f: F)
where
    F: FnMut(&[u8]),
{
    let n = s.len();
    if n == 0 {
        return;
    }

    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(s);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| s[(start + k) % n]).collect();

    // suffix array construction
    let mut sa = vec![0; aperiodic_least_rotation.len()];

    saca::optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);

    println!("Suffix Array: {:?}", sa);

    // rank array
    let rank = rank::sa_to_rank(&sa);
    println!("Rank Array: {:?}", rank);

    // lyndon table
    let lyndon_table = lyndon::lyndontables::lyndon_table(&aperiodic_least_rotation, &rank);
    println!("Lyndon Table: {:?}", lyndon_table);

    // lyndon suffix table
    let lyndon_suffix_table = lyndon::lyndontables::lyndon_suffix_table(&aperiodic_least_rotation);
    println!("Lyndon Suffix Table: {:?}", lyndon_suffix_table);

    for i in 1..period_length {
        let factor_start = i;
        let factor_length = lyndon_table[i];
        let factor = &aperiodic_least_rotation[factor_start..factor_start + factor_length];
        f(factor);
    }

    for i in 0..period_length-1 {
        let factor_end = i + 1;
        let factor_length = lyndon_suffix_table[i];
        let factor_start = factor_end - factor_length;
        let factor = &aperiodic_least_rotation[factor_start..factor_end];
        f(factor);
    }

}

// returns each factor of each conjugate of the string,
//  together with a list of counts of how many times the factor appears in each conjugate
pub fn get_factor_with_conjugates_counts(aperiodic_least_rotation: &[u8]) -> std::collections::HashMap<Vec<u8>, Vec<usize>> {
    let n = aperiodic_least_rotation.len();
    if n == 0 {
        return std::collections::HashMap::new();
    }

    // suffix array construction
    let mut sa = vec![0; n];
    saca::optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);
    println!("Suffix Array: {:?}", sa);

    // rank array
    let rank = rank::sa_to_rank(&sa);
    println!("Rank Array: {:?}", rank);

    // lyndon table
    let lyndon_table = lyndon::lyndontables::lyndon_table(&aperiodic_least_rotation, &rank);
    println!("Lyndon Table: {:?}", lyndon_table);

    // lyndon suffix table
    let lyndon_suffix_table = lyndon::lyndontables::lyndon_suffix_table(&aperiodic_least_rotation);
    println!("Lyndon Suffix Table: {:?}", lyndon_suffix_table);

    let mut factors_with_counts: std::collections::HashMap<Vec<u8>, Vec<usize>> = std::collections::HashMap::new();
    
    // the whole string is a factor that appears only in the first conjugate
    factors_with_counts.insert(aperiodic_least_rotation.to_vec(), vec![0; n]); 
    factors_with_counts.get_mut(aperiodic_least_rotation).unwrap()[0] = 1;

    // Process suffixes
    for i in (1..n).rev() {
        let mut factor_start = i;
        
        while factor_start < n {
            let factor_length = lyndon_table[factor_start];
            let factor = &aperiodic_least_rotation[factor_start..factor_start + factor_length];
            
            // if the factor is missing, add the couple (factor, empty vector of n zeros) to the set
            if !factors_with_counts.contains_key(factor) {
                factors_with_counts.insert(factor.to_vec(), vec![0; n]);
            }

            // Increment the count for this factor in the conjugate starting at index i
            factors_with_counts.get_mut(factor).unwrap()[i] += 1;

            factor_start += factor_length;
        }
    }

    // Process prefixes
    for i in 1..n {
        let mut my_factor_end = i;
        
        while my_factor_end > 0 {
            // FIX: Use lyndon_suffix_table instead of lyndon_table
            let my_factor_length = lyndon_suffix_table[my_factor_end - 1];
            let factor_start = my_factor_end - my_factor_length;
            let factor = &aperiodic_least_rotation[factor_start..my_factor_end];
            
            // if the factor is missing, add the couple (factor, empty vector of n zeros) to the set
            if !factors_with_counts.contains_key(factor) {
                factors_with_counts.insert(factor.to_vec(), vec![0; n]);
            }

            // Increment the count for this factor in the conjugate starting at index i
            factors_with_counts.get_mut(factor).unwrap()[i] += 1;

            // Safely step backward to prevent usize underflow
            my_factor_end = factor_start;
        }
    }

    factors_with_counts
}

pub fn print_all_factors(s: &[u8]) {
    let n = s.len();
    if n == 0 {
        println!("Empty string: nothing to factor.");
        return;
    }

    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(s);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| s[(start + k) % n]).collect();

    println!("Aperiodic least rotation: {}", String::from_utf8_lossy(&aperiodic_least_rotation));
    println!("Period length: {}", period_length);
    println!(
        "So each letter in the BBWT of the string is repeated {} times.",
        n / period_length
    );


    // suffix array construction
    let mut sa = vec![0; aperiodic_least_rotation.len()];

    saca::optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);

    println!("Suffix Array: {:?}", sa);

    // rank array
    let rank = rank::sa_to_rank(&sa);
    println!("Rank Array: {:?}", rank);

    // lyndon table
    let lyndon_table = lyndon::lyndontables::lyndon_table(&aperiodic_least_rotation, &rank);
    println!("Lyndon Table: {:?}", lyndon_table);

    // lyndon suffix table
    let lyndon_suffix_table = lyndon::lyndontables::lyndon_suffix_table(&aperiodic_least_rotation);
    println!("Lyndon Suffix Table: {:?}", lyndon_suffix_table);

    // heuristics
    let most_number_of_factors_index = crate::heuristics::most_factors(&lyndon_table, &lyndon_suffix_table);
    println!("Maximum number of factors: {}", most_number_of_factors_index);

    let most_repeated_factors_index = crate::heuristics::most_repeated_factors(
        &aperiodic_least_rotation,
        &lyndon_table,
        &lyndon_suffix_table,
    );
    println!("Index with most repeated factors: {}", most_repeated_factors_index);


    // print BBWT for each heuristic and the original string
    let bwt = crate::bbwt::bbwt(&aperiodic_least_rotation);
    println!("BWT: {}", String::from_utf8_lossy(&bwt));

    let bwt_most_factors = crate::bbwt::bbwt(&rotate_string(s, most_number_of_factors_index));
    println!("BWT (most factors): {}", String::from_utf8_lossy(&bwt_most_factors));

    let bwt_most_repeated_factors = crate::bbwt::bbwt(&rotate_string(s, most_repeated_factors_index));
    println!("BWT (most repeated factors): {}", String::from_utf8_lossy(&bwt_most_repeated_factors));

    // lcp
    let lcp_array = lyndon::lcp::lcp(&aperiodic_least_rotation, &sa, &rank);
    println!("LCP Array: {:?}", lcp_array);

    online_bbwt_check(&aperiodic_least_rotation);

}