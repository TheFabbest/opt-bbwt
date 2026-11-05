use crate::lyndon::{self, shiloach};
use crate::suffix_array::suffix_array;
use crate::{rank};
use rustc_hash::FxHashMap;

// takes a lambda and executes it on each factor of each conjugate of the string
pub fn for_each_factor_of_each_conjugate<F>(s: &[u8], mut f: F)
where
    F: FnMut(&[u8]),
{
    let n = s.len();
    if n == 0 {
        return;
    }

    let (start, period_length) = shiloach(s);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| s[(start + k) % n]).collect();

    // suffix array construction
    let sa = suffix_array(&aperiodic_least_rotation);

    // rank array
    let rank = rank::sa_to_rank(&sa);

    // lyndon table
    let lyndon_table = lyndon::lyndontables::lyndon_table(&aperiodic_least_rotation, &rank);

    // lyndon suffix table
    let lyndon_suffix_table = lyndon::lyndontables::lyndon_suffix_table(&aperiodic_least_rotation);

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


// A completely flat structure designed for cache-friendly iteration
pub struct GroupedFactors<'a> {
    pub sorted_factors: Vec<&'a [u8]>,
    pub factor_indices: Vec<usize>,       // Flat array of all conjugate indices
    pub factor_starts: Vec<usize>,        // Start offset in factor_indices for each factor
    pub factor_lengths: Vec<usize>,       // Number of indices for each factor
}

pub fn get_factor_with_conjugates_grouped<'a>(
    aperiodic_least_rotation: &'a [u8],
) -> GroupedFactors<'a> {
    let n = aperiodic_least_rotation.len();
    if n == 0 {
        return GroupedFactors {
            sorted_factors: Vec::new(),
            factor_indices: Vec::new(),
            factor_starts: Vec::new(),
            factor_lengths: Vec::new(),
        };
    }

    let sa = suffix_array(aperiodic_least_rotation);

    let rank = rank::sa_to_rank(&sa);
    let lyndon_table = lyndon::lyndontables::lyndon_table(aperiodic_least_rotation, &rank);
    let lyndon_suffix_table = lyndon::lyndontables::lyndon_suffix_table(aperiodic_least_rotation);

    // Arena for linked list: stores (conjugate_index, next_node_index)
    let mut arena: Vec<(usize, usize)> = Vec::with_capacity(n * 2);
    // Map: factor -> (head_idx, tail_idx)
    let mut factors_map: FxHashMap<&'a [u8], (usize, usize)> = FxHashMap::default();

    let mut push_to_arena = |factor: &'a [u8], idx: usize| {
        use std::collections::hash_map::Entry;
        match factors_map.entry(factor) {
            Entry::Occupied(mut occ) => {
                let tail = &mut occ.get_mut().1;
                // Prevent basic consecutive duplicates
                if arena[*tail].0 != idx {
                    let new_node = arena.len();
                    arena.push((idx, usize::MAX));
                    arena[*tail].1 = new_node;
                    *tail = new_node; // Update tail in place
                }
            }
            Entry::Vacant(vac) => {
                let new_node = arena.len();
                arena.push((idx, usize::MAX));
                vac.insert((new_node, new_node));
            }
        }
    };

    push_to_arena(aperiodic_least_rotation, 0);

    // Process suffixes
    for i in (1..n).rev() {
        let mut factor_start = i;
        while factor_start < n {
            let factor_length = lyndon_table[factor_start];
            let factor = &aperiodic_least_rotation[factor_start..factor_start + factor_length];
            push_to_arena(factor, i);
            factor_start += factor_length;
        }
    }

    // Process prefixes
    for i in 1..n {
        let mut my_factor_end = i;
        while my_factor_end > 0 {
            let my_factor_length = lyndon_suffix_table[my_factor_end - 1];
            let factor_start = my_factor_end - my_factor_length;
            let factor = &aperiodic_least_rotation[factor_start..my_factor_end];
            push_to_arena(factor, i);
            my_factor_end = factor_start;
        }
    }

    // 1. Collect unique factors and sort them reverse-lexicographically
    let mut sorted_factors: Vec<&'a [u8]> = factors_map.keys().copied().collect();
    sorted_factors.sort_unstable_by(|a, b| b.cmp(a));

    let mut factor_indices = Vec::with_capacity(arena.len());
    let mut factor_starts = Vec::with_capacity(sorted_factors.len());
    let mut factor_lengths = Vec::with_capacity(sorted_factors.len());

    // 2. Flatten the linked lists into a single continuous array
    for factor in &sorted_factors {
        let (head, _) = factors_map[factor];
        let start = factor_indices.len();
        
        let mut curr = head;
        while curr != usize::MAX {
            factor_indices.push(arena[curr].0);
            curr = arena[curr].1;
        }
        
        // Sort the newly added slice
        factor_indices[start..].sort_unstable();
        
        // Stable in-place deduplication for the tail of the vector
        let mut write_idx = start;
        for read_idx in start..factor_indices.len() {
            if write_idx == start || factor_indices[read_idx] != factor_indices[write_idx - 1] {
                factor_indices[write_idx] = factor_indices[read_idx];
                write_idx += 1;
            }
        }
        
        // Truncate back to the deduplicated length
        factor_indices.truncate(write_idx);
        let keep_len = write_idx - start;

        factor_starts.push(start);
        factor_lengths.push(keep_len);
    }

    GroupedFactors {
        sorted_factors,
        factor_indices,
        factor_starts,
        factor_lengths,
    }
}