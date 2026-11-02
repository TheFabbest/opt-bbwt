const L_TYPE: bool = false;
const S_TYPE: bool = true;

#[inline]
fn is_lms(types: &[bool], i: usize) -> bool {
    i > 0 && types[i] == S_TYPE && types[i - 1] == L_TYPE
}

fn type_array(s: &[usize]) -> Vec<bool> {
    let mut types = vec![L_TYPE; s.len()];
    types[s.len() - 1] = S_TYPE; // Sentinels are always S-type
    for i in (0..s.len() - 1).rev() {
        if s[i] < s[i + 1] {
            types[i] = S_TYPE;
        } else if s[i] == s[i + 1] {
            types[i] = types[i + 1];
        } else {
            types[i] = L_TYPE;
        }
    }
    types
}

fn bucket_counts(s: &[usize], alphabet_size: usize) -> Vec<usize> {
    let mut counts = vec![0; alphabet_size];
    for &c in s {
        counts[c] += 1;
    }
    counts
}

fn bucket_heads(counts: &[usize]) -> Vec<usize> {
    let mut heads = vec![0; counts.len()];
    let mut sum = 0;
    for i in 0..counts.len() {
        heads[i] = sum;
        sum += counts[i];
    }
    heads
}

fn bucket_tails(counts: &[usize]) -> Vec<usize> {
    let mut tails = vec![0; counts.len()];
    let mut sum = 0;
    for i in 0..counts.len() {
        sum += counts[i];
        tails[i] = sum;
    }
    tails
}

fn induce_sort_l(s: &[usize], sa: &mut [usize], heads: &[usize], types: &[bool]) {
    let mut head_ptrs = heads.to_vec();
    for i in 0..sa.len() {
        if sa[i] != usize::MAX && sa[i] > 0 {
            let j = sa[i] - 1;
            if types[j] == L_TYPE {
                let c = s[j];
                sa[head_ptrs[c]] = j;
                head_ptrs[c] += 1;
            }
        }
    }
}

fn induce_sort_s(s: &[usize], sa: &mut [usize], tails: &[usize], types: &[bool]) {
    let mut tail_ptrs = tails.to_vec();
    for i in (0..sa.len()).rev() {
        if sa[i] != usize::MAX && sa[i] > 0 {
            let j = sa[i] - 1;
            if types[j] == S_TYPE {
                let c = s[j];
                tail_ptrs[c] -= 1;
                sa[tail_ptrs[c]] = j;
            }
        }
    }
}

fn are_lms_substrings_equal(s: &[usize], types: &[bool], i: usize, j: usize) -> bool {
    if i == s.len() || j == s.len() {
        return i == j;
    }
    let mut k = 0;
    loop {
        let i_lms = is_lms(types, i + k);
        let j_lms = is_lms(types, j + k);

        if k > 0 && i_lms && j_lms {
            return true;
        }
        if i_lms != j_lms || s[i + k] != s[j + k] {
            return false;
        }
        k += 1;
        if i + k >= s.len() || j + k >= s.len() {
            return false;
        }
    }
}

fn sais_core(s: &[usize], alphabet_size: usize) -> Vec<usize> {
    let len = s.len();
    if len == 1 {
        return vec![0];
    }

    let types = type_array(s);
    let counts = bucket_counts(s, alphabet_size);
    let heads = bucket_heads(&counts);
    let tails = bucket_tails(&counts);

    let mut sa = vec![usize::MAX; len];

    // Step 1: Guess LMS positions
    let mut tail_ptrs = tails.clone();
    for i in (0..len).rev() {
        if is_lms(&types, i) {
            let c = s[i];
            tail_ptrs[c] -= 1;
            sa[tail_ptrs[c]] = i;
        }
    }

    // Step 2: Induce sort
    induce_sort_l(s, &mut sa, &heads, &types);
    induce_sort_s(s, &mut sa, &tails, &types);

    // Step 3: Summarize LMS substrings
    let mut lms_map = vec![usize::MAX; len];
    let mut lms_offsets = Vec::new();
    
    // SAFE LMS EXTRACTION: Collect LMS offsets in original string order
    for i in 0..len {
        if is_lms(&types, i) {
            lms_map[i] = lms_offsets.len();
            lms_offsets.push(i);
        }
    }

    let lms_count = lms_offsets.len();
    if lms_count == 0 {
        return sa;
    }

    let mut summary_s = vec![0; lms_count];
    let mut current_name = 0;
    let mut prev_lms = usize::MAX;

    for &i in sa.iter() {
        if i != usize::MAX && is_lms(&types, i) {
            if prev_lms != usize::MAX && !are_lms_substrings_equal(s, &types, prev_lms, i) {
                current_name += 1;
            }
            summary_s[lms_map[i]] = current_name;
            prev_lms = i;
        }
    }

    // Step 4: Solve summary SA
    let summary_sa = if current_name + 1 == lms_count {
        let mut sa_sum = vec![0; lms_count];
        for (i, &name) in summary_s.iter().enumerate() {
            sa_sum[name] = i;
        }
        sa_sum
    } else {
        sais_core(&summary_s, current_name + 1)
    };

    // Step 5: Accurate LMS sorting & final induction
    sa.fill(usize::MAX);
    let mut tail_ptrs = tails.clone();
    for &idx in summary_sa.iter().rev() {
        let string_index = lms_offsets[idx];
        let c = s[string_index];
        tail_ptrs[c] -= 1;
        sa[tail_ptrs[c]] = string_index;
    }

    induce_sort_l(s, &mut sa, &heads, &types);
    induce_sort_s(s, &mut sa, &tails, &types);

    sa
}

pub fn sais(s: &[u8]) -> Vec<usize> {
    if s.is_empty() {
        return Vec::new();
    }

    let mut max_val = 0;
    for &c in s {
        if (c as usize) > max_val {
            max_val = c as usize;
        }
    }

    // Map input to internal alphabet reserving '0' as sentinel '$'
    let mut mapped = Vec::with_capacity(s.len() + 1);
    for &c in s {
        mapped.push((c as usize) + 1);
    }
    mapped.push(0); // Virtual sentinel '$'

    let full_sa = sais_core(&mapped, max_val + 2);

    // Omit the sentinel suffix from output
    full_sa[1..].to_vec()
}