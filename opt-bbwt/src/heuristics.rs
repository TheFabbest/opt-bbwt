//! This module contains functions that select the conjugate of a string that shows a certain property that usually leads to a BBWT with fewer runs.
//! These heuristics are not guaranteed to find the optimal conjugate, but they are fast and often find a good one.
//! More details can be found in <a href="https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates">Bijective BWT: Implementing Linear-Time Lyndon Factorization of Conjugates and Discussing Open Problems</a>.

pub fn most_factors(lyndon_table: &[usize], lyndon_suffix_table: &[usize]) -> usize {
    let len = lyndon_table.len();
    if len == 0 {
        return 0;
    }

    let mut counts_for_suffixes = vec![0; len + 1];
    counts_for_suffixes[0] = 1;
    counts_for_suffixes[len] = 0;
    for i in (1..len).rev() {
        counts_for_suffixes[i] = 1 + counts_for_suffixes[i + lyndon_table[i]];
    }

    let mut counts_for_prefixes = vec![0; len + 1];
    counts_for_prefixes[0] = 0;
    counts_for_prefixes[len] = 1;
    for i in 0..len - 1 {
        counts_for_prefixes[i + 1] = 1 + counts_for_prefixes[i + 1 - lyndon_suffix_table[i]];
    }

    let mut max_factors = 0;
    let mut max_index = 0;
    for i in 1..len {
        let total_factors = counts_for_suffixes[i] + counts_for_prefixes[i];
        if total_factors > max_factors {
            max_factors = total_factors;
            max_index = i;
        }
    }
    max_index
}

pub fn most_repeated_factors(
    s: &[u8],
    lyndon_table: &[usize],
    lyndon_suffix_table: &[usize],
) -> usize {
    compute_repeated_factors(s, lyndon_table, lyndon_suffix_table, |_| 1)
}

pub fn most_repeated_factors_proportional(
    s: &[u8],
    lyndon_table: &[usize],
    lyndon_suffix_table: &[usize],
) -> usize {
    compute_repeated_factors(s, lyndon_table, lyndon_suffix_table, |length| length)
}

fn compute_repeated_factors(
    s: &[u8],
    lyndon_table: &[usize],
    lyndon_suffix_table: &[usize],
    weight: impl Fn(usize) -> usize,
) -> usize {
    let len = lyndon_table.len();

    let mut repetitions_in_suffixes = vec![0; len + 1];
    repetitions_in_suffixes[len] = 0;
    for i in (1..len).rev() {
        let my_factor_start = i;
        let my_factor_end = i + lyndon_table[i];
        let next_factor_start = my_factor_end;
        if next_factor_start >= len {
            repetitions_in_suffixes[i] = 0;
            continue;
        }
        let next_factor_end = next_factor_start + lyndon_table[next_factor_start];

        if s[my_factor_start..my_factor_end] == s[next_factor_start..next_factor_end] {
            repetitions_in_suffixes[i] = 1 + repetitions_in_suffixes[i + lyndon_table[i]];
        } else {
            repetitions_in_suffixes[i] = repetitions_in_suffixes[i + lyndon_table[i]];
        }
    }

    let mut repetitions_in_prefixes = vec![0; len + 1];
    repetitions_in_prefixes[0] = 0;
    for i in 0..len - 1 {
        let my_factor_start = i + 1 - lyndon_suffix_table[i];
        if my_factor_start == 0 {
            repetitions_in_prefixes[i + 1] = 0;
            continue;
        }
        let my_factor_end = i + 1;
        let next_factor_start = my_factor_start - lyndon_suffix_table[my_factor_start - 1];
        let next_factor_end = my_factor_start;

        if s[my_factor_start..my_factor_end] == s[next_factor_start..next_factor_end] {
            repetitions_in_prefixes[i + 1] =
                weight(my_factor_end - my_factor_start) + repetitions_in_prefixes[my_factor_start];
        } else {
            repetitions_in_prefixes[i + 1] = repetitions_in_prefixes[my_factor_start];
        }
    }

    let mut max_repeated_factors = 0;
    let mut max_index = 0;
    for i in 1..len {
        let total_repetitions = repetitions_in_suffixes[i] + repetitions_in_prefixes[i];
        if total_repetitions > max_repeated_factors {
            max_repeated_factors = total_repetitions;
            max_index = i;
        }
    }
    max_index
}

#[cfg(test)]
mod tests {

    use crate::{
        file::generate_random_binary_string,
        lyndon::{duval, lyndon_suffix_table, lyndon_table, shiloach},
        string::rotate_string,
        suffix_array::{sa_to_rank, suffix_array},
    };

    use super::*;

    #[test]
    fn test_most_factors_empty() {
        let lyndon_table = vec![];
        let lyndon_suffix_table = vec![];
        let result = most_factors(&lyndon_table, &lyndon_suffix_table);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_most_factors_single_char() {
        let lyndon_table = vec![1];
        let lyndon_suffix_table = vec![1];
        let result = most_factors(&lyndon_table, &lyndon_suffix_table);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_most_factors_example() {
        let lyndon_table = vec![6, 1, 2, 1, 2, 1];
        let lyndon_suffix_table = vec![1, 2, 1, 4, 1, 6];
        let result = most_factors(&lyndon_table, &lyndon_suffix_table);
        assert!(result == 1 || result == 3);
    }

    #[test]
    fn test_most_repeated_factors_example() {
        let s = b"abanan";
        let lyndon_table = vec![6, 1, 2, 1, 2, 1];
        let lyndon_suffix_table = vec![1, 2, 1, 4, 1, 6];
        let result = most_repeated_factors(s, &lyndon_table, &lyndon_suffix_table);
        assert!(result == 1 || result == 2);
    }

    #[test]
    fn test_most_repeated_factors_random() {
        let n = 100;
        let text = generate_random_binary_string(n).into_bytes();
        let (start, period_length) = shiloach(&text);
        let s: Vec<u8> = (0..period_length).map(|k| text[(start + k) % n]).collect();

        let sa = suffix_array(&s);
        let rank = sa_to_rank(&sa);
        let lyndon_table = lyndon_table(&s, &rank);
        let lyndon_suffix_table = lyndon_suffix_table(&s);
        let most: usize = most_repeated_factors(&s, &lyndon_table, &lyndon_suffix_table);
        let found_repetitions = duval(&rotate_string(&s, most))
            .windows(2)
            .filter(|w| w[0] == w[1])
            .count();

        for i in 0..s.len() - 1 {
            let rotated = rotate_string(&s, i);
            let factors = duval(&rotated);
            let repetitions = factors.windows(2).filter(|w| w[0] == w[1]).count();
            assert!(
                repetitions <= found_repetitions,
                "Repetitions: {}, Max Repetitions: {}, Index: {}",
                repetitions,
                found_repetitions,
                i
            );
        }
    }

    #[test]
    fn test_most_repeated_factors_edge_characters() {
        let s = (0u8..=255).collect::<Vec<u8>>();
        let (start, period_length) = shiloach(&s);
        let s: Vec<u8> = (0..period_length)
            .map(|k| s[(start + k) % s.len()])
            .collect();

        let sa = suffix_array(&s);
        let rank = sa_to_rank(&sa);
        let lyndon_table = lyndon_table(&s, &rank);
        let lyndon_suffix_table = lyndon_suffix_table(&s);
        let most: usize = most_repeated_factors(&s, &lyndon_table, &lyndon_suffix_table);
        let found_repetitions = duval(&rotate_string(&s, most))
            .windows(2)
            .filter(|w| w[0] == w[1])
            .count();

        for i in 0..s.len() - 1 {
            let rotated = rotate_string(&s, i);
            let factors = duval(&rotated);
            let repetitions = factors.windows(2).filter(|w| w[0] == w[1]).count();
            assert!(
                repetitions <= found_repetitions,
                "Repetitions: {}, Max Repetitions: {}, Index: {}",
                repetitions,
                found_repetitions,
                i
            );
        }
    }
}
