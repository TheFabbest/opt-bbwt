use crate::{lyndon::{duval, lyndon_suffix_table, lyndon_table, shiloach}, rank::sa_to_rank, suffix_array::suffix_array, string::from_period_slice_into};

const L_TYPE: bool = false;
const S_TYPE: bool = true;

/// Helper BitVector-like wrapper to mirror C++ BitVector functionality
struct BitVector {
    bits: Vec<bool>,
}

impl BitVector {
    fn new(len: usize) -> Self {
        Self {
            bits: vec![false; len],
        }
    }

    fn get(&self, pos: usize) -> bool {
        if pos < self.bits.len() {
            self.bits[pos]
        } else {
            false
        }
    }

    fn set(&mut self, pos: usize, val: bool) {
        if pos < self.bits.len() {
            self.bits[pos] = val;
        }
    }

    /// Finds next position STRICTLY GREATER THAN pos where bit is true
    fn next(&self, pos: usize) -> usize {
        for i in (pos + 1)..self.bits.len() {
            if self.bits[i] {
                return i;
            }
        }
        self.bits.len()
    }

    /// Finds previous position STRICTLY LESS THAN pos where bit is true
    fn prev(&self, pos: isize) -> isize {
        let mut p = pos - 1;
        let max_pos = (self.bits.len() - 1) as isize;
        if p > max_pos {
            p = max_pos;
        }
        
        while p >= 0 {
            if self.bits[p as usize] {
                return p;
            }
            p -= 1;
        }
        -1
    }
}

fn is_lms_pos(pos: usize, l_fac: &BitVector, s_type: &BitVector) -> bool {
    if pos == 0 {
        return l_fac.get(0);
    }
    l_fac.get(pos) || (s_type.get(pos) == S_TYPE && s_type.get(pos - 1) == L_TYPE)
}

fn compute_buckets_structure(in_str: &[usize], buckets: &mut [usize], alph_size: usize) {
    let mut char_count = vec![0; alph_size];
    for &c in in_str {
        char_count[c] += 1;
    }

    buckets[0] = 0;
    let mut total = 0;
    for i in 0..alph_size {
        total += char_count[i];
        buckets[i + 1] = total;
    }
}

fn pre_sort_suffixes_l(
    in_str: &[usize],
    sa: &mut [isize],
    len: usize,
    l_fac: &BitVector,
    suff_type: &BitVector,
    spc_fac: &BitVector,
    buckets: &mut [usize],
) {
    let mut p = spc_fac.prev(len as isize);

    for i in 0..len {
        while p >= 0 && buckets[in_str[p as usize]] == i {
            let j = l_fac.next(p as usize) - 1;
            sa[buckets[in_str[j]]] = j as isize;
            buckets[in_str[j]] += 1;
            p = spc_fac.prev(p);
        }

        let j = sa[i];
        if j < 0 {
            continue;
        }
        let mut j = j as usize;

        if !l_fac.get(j) {
            j = j - 1;
        } else {
            j = l_fac.next(j) - 1;
        }

        if suff_type.get(j) == L_TYPE {
            sa[buckets[in_str[j]]] = j as isize;
            buckets[in_str[j]] += 1;
        }
    }
}

fn pre_sort_suffixes_s(
    in_str: &[usize],
    sa: &mut [isize],
    len: usize,
    l_fac: &BitVector,
    suff_type: &BitVector,
    buckets: &mut [usize],
) {
    for i in (0..len).rev() {
        let j = sa[i];
        if j < 0 {
            continue;
        }

        let mut j = j as usize;
        if !l_fac.get(j) {
            j = j - 1;

            if suff_type.get(j) == S_TYPE {
                buckets[in_str[j] + 1] -= 1;
                sa[buckets[in_str[j] + 1]] = j as isize;
            }
        }
    }
}

pub fn circular_suffix_array(text: &[u8]) -> Vec<usize> {
    let mut output = vec![0; text.len()];
    circular_suffix_array_into(text, &mut output);
    output
}

/// Computes the circular suffix array into a caller-provided buffer.
///
/// `output` must have the same length as `text`.
///
/// # Panics
///
/// Panics if `output.len() != text.len()`.
pub fn circular_suffix_array_into(text: &[u8], output: &mut [usize]) {
    assert_eq!(
        output.len(),
        text.len(),
        "circular suffix-array output buffer must match input length"
    );
    let n = text.len();
    if n == 0 {
        return;
    }
    let min_c = *text.iter().min().unwrap() as usize;
    let max_c = *text.iter().max().unwrap() as usize;
    let alph_size = max_c - min_c + 1;

    let in_str: Vec<usize> = text.iter().map(|&c| (c as usize) - min_c).collect();
    let factors = duval(text);

    let mut lb_fac = BitVector::new(n + 1);
    for &f in &factors {
        lb_fac.set(f, true);
    }
    lb_fac.set(n, true);

    let mut buckets = vec![0; alph_size + 1];
    compute_buckets_structure(&in_str, &mut buckets, alph_size);

    let mut sa = vec![-1isize; n];
    circular_suffix_array_inner(&in_str, &mut sa, n, &lb_fac, alph_size, &buckets);

    for (output_index, value) in sa.into_iter().enumerate() {
        output[output_index] = value as usize;
    }
}

fn circular_suffix_array_inner(
    in_str: &[usize],
    sa: &mut [isize],
    len: usize,
    lb_fac: &BitVector,
    alph_size: usize,
    buckets: &Vec<usize>,
) {
    if len == 0 {
        return;
    }

    let mut suff_type = BitVector::new(len + 7);
    let mut spc_suff = BitVector::new(len + 1);

    let mut f_start = 0;
    while f_start < len {
        let f_end = lb_fac.next(f_start);
        suff_type.set(f_start, S_TYPE);

        if f_end >= 2 {
            for j in (f_start..=f_end - 2).rev() {
                if (in_str[j] < in_str[j + 1])
                    || (in_str[j] == in_str[j + 1] && suff_type.get(j + 1) == S_TYPE)
                {
                    suff_type.set(j, S_TYPE);
                }
            }
        }

        let mut m = 0;
        let mut c = 0;
        let mut c1 = in_str[f_end - 1];

        if f_end >= 2 {
            for i in (f_start..=f_end - 2).rev() {
                let c0 = in_str[i];
                if c0 < c1 + c {
                    c = 1;
                } else if c != 0 {
                    m += 1;
                    c = 0;
                }
                c1 = c0;
            }
        }

        if m == 0 && c == 0 {
            spc_suff.set(f_start, true);
        }

        f_start = f_end;
    }

    spc_suff.set(len, true);

    
    let mut tmp_buckets = vec![0; alph_size + 1];

    for i in 0..len {
        sa[i] = -1;
    }

    tmp_buckets.copy_from_slice(&buckets);
    for i in 0..len {
        if is_lms_pos(i, lb_fac, &suff_type) && !spc_suff.get(i) {
            let b_idx = tmp_buckets[in_str[i] + 1] - 1;
            sa[b_idx] = i as isize;
            tmp_buckets[in_str[i] + 1] -= 1;
        }
    }

    tmp_buckets.copy_from_slice(&buckets);
    pre_sort_suffixes_l(in_str, sa, len, lb_fac, &suff_type, &spc_suff, &mut tmp_buckets);

    tmp_buckets.copy_from_slice(&buckets);
    pre_sort_suffixes_s(in_str, sa, len, lb_fac, &suff_type, &mut tmp_buckets);

    let mut num_lms_suff = 0;
    for i in 0..len {
        let pos = sa[i];
        if pos >= 0 && is_lms_pos(pos as usize, lb_fac, &suff_type) && !spc_suff.get(pos as usize) {
            sa[num_lms_suff] = pos;
            num_lms_suff += 1;
        }
    }

    for i in num_lms_suff..len {
        sa[i] = 0;
    }

    f_start = 0;
    while f_start < len {
        let f_end = lb_fac.next(f_start);
        let mut j = f_end;
        let mut c = 0;
        let mut c1 = in_str[f_end - 1];

        if f_end >= 2 {
            for i in (f_start..=f_end - 2).rev() {
                let c0 = in_str[i];
                if c0 < c1 + c {
                    c = 1;
                } else if c != 0 {
                    sa[num_lms_suff + ((i + 1) >> 1)] = (j - i - 1) as isize;
                    j = i + 1;
                    c = 0;
                }
                c1 = c0;
            }
        }

        if j < f_end || c != 0 {
            sa[num_lms_suff + (f_start >> 1)] = (j - f_start) as isize;
        }

        f_start = f_end;
    }

    let mut num_labels = 0;
    let mut q = len;
    let mut q_len = 0;

    for i in 0..num_lms_suff {
        let pos = sa[i] as usize;
        let sbwrd_len = sa[num_lms_suff + (pos >> 1)] as usize;
        let mut distinct = true;

        if sbwrd_len == q_len {
            let mut j = 0;
            while j < sbwrd_len && in_str[pos + j] == in_str[q + j] {
                j += 1;
            }
            if j == sbwrd_len {
                distinct = false;
            }
        }

        if distinct {
            num_labels += 1;
            q = pos;
            q_len = sbwrd_len;
        }

        sa[num_lms_suff + (pos >> 1)] = num_labels as isize;
    }

    if num_labels < num_lms_suff {
        let mut red_factors = BitVector::new(num_lms_suff + 1);

        let mut out_pos = 0;
        for in_pos in 0..len {
            if is_lms_pos(in_pos, lb_fac, &suff_type) && !spc_suff.get(in_pos) {
                red_factors.set(out_pos, lb_fac.get(in_pos));
                out_pos += 1;
            }
        }
        red_factors.set(num_lms_suff, true);

        let mut red_str = vec![0usize; num_lms_suff];
        let mut out_pos = num_lms_suff;
        for in_pos in (num_lms_suff..len).rev() {
            if sa[in_pos] != 0 {
                out_pos -= 1;
                red_str[out_pos] = (sa[in_pos] - 1) as usize;
            }
        }

        let mut rec_sa = vec![-1isize; num_lms_suff];
        let mut new_buckets = vec![0; num_labels + 1];
        compute_buckets_structure(&red_str, &mut new_buckets, num_labels);
        circular_suffix_array_inner(&red_str, &mut rec_sa, num_lms_suff, &red_factors, num_labels, &new_buckets);

        out_pos = 0;
        for in_pos in 0..len {
            if is_lms_pos(in_pos, lb_fac, &suff_type) && !spc_suff.get(in_pos) {
                red_str[out_pos] = in_pos;
                out_pos += 1;
            }
        }

        for i in 0..num_lms_suff {
            sa[i] = red_str[rec_sa[i] as usize] as isize;
        }
    }

    tmp_buckets.copy_from_slice(&buckets);
    for i in num_lms_suff..len {
        sa[i] = -1;
    }

    for i in (0..num_lms_suff).rev() {
        let j = sa[i] as usize;
        sa[i] = -1;
        let b_idx = tmp_buckets[in_str[j] + 1] - 1;
        sa[b_idx] = j as isize;
        tmp_buckets[in_str[j] + 1] -= 1;
    }

    tmp_buckets.copy_from_slice(&buckets);
    pre_sort_suffixes_l(in_str, sa, len, lb_fac, &suff_type, &spc_suff, &mut tmp_buckets);

    tmp_buckets.copy_from_slice(&buckets);
    pre_sort_suffixes_s(in_str, sa, len, lb_fac, &suff_type, &mut tmp_buckets);
}


pub fn optimal_csa(input: &[u8]) -> usize {
    let (start, period_length) = shiloach(input);

    let mut text_storage = vec![0; period_length];
    from_period_slice_into(input, start, period_length, &mut text_storage);
    let text = text_storage.as_slice();


    let min_c = *text.iter().min().unwrap() as usize;
    let max_c = *text.iter().max().unwrap() as usize;
    let alph_size = max_c - min_c + 1;

    let mut buckets = vec![0; alph_size + 1];

    let in_str: Vec<usize> = text.iter().map(|&c| (c as usize) - min_c).collect();

    compute_buckets_structure(&in_str, &mut buckets, alph_size);

    let n = in_str.len();

    let mut concatenated = Vec::with_capacity(2 * n);
    concatenated.extend_from_slice(&in_str);
    concatenated.extend_from_slice(&in_str);

    let sa = &suffix_array(text);
    let rank = &sa_to_rank(sa);

    let lyndon_table = lyndon_table(text, rank);
    let lyndon_suffix_table = lyndon_suffix_table(text);

    let mut best_run_count = 0;
    let mut best_conjugate = 0;
    let mut csa = vec![0; n];

    let mut prev_char = None;
    for idx in sa {
        let prev_idx = if *idx == 0 {
            n - 1
        } else {
            *idx - 1
        };

        let current = text[prev_idx];

        if prev_char != Some(current) {
            best_run_count += 1;
            prev_char = Some(current);
        }
    }

    for shift in 1..n {
        let start = (n - shift) % n;
        let rotated_string = &concatenated[n - start.. n + n - start].to_vec();

        
        let mut factors: Vec<usize> = Vec::new();

        let mut factor_start = shift;
        while factor_start < n {
            let factor_len = lyndon_table[factor_start];
            factors.push(factor_start - shift);
            factor_start += factor_len;
        }
        
        let mut u_factor_lengths = Vec::new();
        let mut end = shift;
        while end > 0 {
            let factor_len = lyndon_suffix_table[end - 1];
            u_factor_lengths.push(factor_len);
            end -= factor_len;
        }

        u_factor_lengths.reverse();

        let mut current_offset = n - shift;
        for len in u_factor_lengths {
            factors.push(current_offset);
            current_offset += len;
        }


        let mut lb_fac = BitVector::new(n + 1);
        for &f in &factors {
            lb_fac.set(f, true);
        }
        lb_fac.set(n, true);

        circular_suffix_array_inner(&rotated_string, &mut csa, n, &lb_fac, alph_size, &buckets);

        // count runs
        let mut run_count = 0;
        let mut prev_char = None;
        for &idx in &csa {
            let factor_start = factors
                .iter()
                .rev()
                .find(|&&start| start <= idx as usize)
                .copied()
                .unwrap_or(0);
            let factor_len = if let Some(&next_start) = factors.iter().find(|&&s| s > factor_start) {
                next_start - factor_start
            } else {
                n - factor_start
            };

            let prev_idx = if idx as usize == factor_start {
                factor_start + factor_len - 1
            } else {
                idx as usize - 1
            };

            let current = rotated_string[prev_idx];

            if prev_char != Some(current) {
                run_count += 1;
                if run_count >= best_run_count {
                    break; 
                }
                prev_char = Some(current);
            }
        }

        if run_count < best_run_count {
            best_run_count = run_count;
            best_conjugate = shift;
        }
    }

    best_conjugate
}