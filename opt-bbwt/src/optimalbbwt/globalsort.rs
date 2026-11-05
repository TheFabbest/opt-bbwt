use crate::lyndon::conjugatesfactors::get_factor_with_conjugates_grouped;
use crate::lyndon::shiloach;
use crate::saca::circularsais::circular_suffix_array;

pub fn optimal_bbwt_global_sort(text: &[u8]) -> usize {
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