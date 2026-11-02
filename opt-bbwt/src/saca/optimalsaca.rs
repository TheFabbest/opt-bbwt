use std::cmp::Ordering;

// Special sentinel constants
const BT: usize = usize::MAX;
const BH: usize = usize::MAX - 1;
const E: usize = usize::MAX - 2;
const R1: usize = usize::MAX - 3;
const R2: usize = usize::MAX - 4;
#[allow(dead_code)]
const MAX_CHAR: usize = usize::MAX - 5;
#[allow(dead_code)]
const MSB: usize = 1usize << (usize::BITS - 1);

// ----- Auxiliary Functions -----

pub fn chars_to_usize_array(input_chars: &[u8]) -> Vec<usize> {
    // .collect() on a slice iterator automatically pre-allocates exact capacity
    input_chars.iter().map(|&c| c as usize).collect()
}

fn get_s_substring_end(input: &[usize], start: usize, length: usize) -> usize {
    if start + 1 >= length || length < 2 {
        return length;
    }
    let mut j = start + 1;
    let limit = length.saturating_sub(1);

    while j < limit && input[j] >= input[j + 1] {
        j += 1;
    }
    if j == limit {
        return length; // sentinel
    }

    let mut k = j;
    while k > start && input[k] == input[j] {
        k -= 1;
    }

    k + 1
}

fn get_l_substring_end(input: &[usize], start: usize, length: usize) -> usize {
    if start + 1 >= length || length < 2 {
        return length;
    }
    let mut j = start + 1;
    let limit = length.saturating_sub(1);

    while j < limit && input[j] <= input[j + 1] {
        j += 1;
    }

    let mut k = j.saturating_sub(1);
    while k > start && input[k] == input[j] {
        k -= 1;
    }

    k + 1
}

fn compare_substrings(
    input: &[usize],
    mut offset_a: usize,
    mut offset_b: usize,
    length: usize,
    using_l_type: bool,
) -> Ordering {
    let (end_a, end_b) = if using_l_type {
        (
            get_l_substring_end(input, offset_a, length),
            get_l_substring_end(input, offset_b, length),
        )
    } else {
        (
            get_s_substring_end(input, offset_a, length),
            get_s_substring_end(input, offset_b, length),
        )
    };

    while offset_a < end_a && offset_b < end_b {
        if input[offset_a] != input[offset_b] {
            return input[offset_a].cmp(&input[offset_b]);
        }
        offset_a += 1;
        offset_b += 1;
    }

    if offset_a == length && offset_b == length {
        Ordering::Equal
    } else if offset_a == length {
        Ordering::Less
    } else if offset_b == length {
        Ordering::Greater
    } else if input[offset_a] != input[offset_b] {
        input[offset_a].cmp(&input[offset_b])
    } else if offset_a == end_a && offset_b == end_b {
        Ordering::Equal
    } else if offset_a == end_a {
        if using_l_type { Ordering::Less } else { Ordering::Greater }
    } else if using_l_type {
        Ordering::Greater
    } else {
        Ordering::Less
    }
}

fn heapify(sa: &mut [usize], n_s: usize, heap_len: usize, mut index: usize) {
    let length = sa.len();
    if length < n_s { return; }
    let (sa_prefix, sa_suffix) = sa.split_at_mut(length - n_s);
    let array = &mut sa_suffix[..n_s];

    loop {
        let mut biggest = index;
        let left = 2 * index + 1;
        let right = 2 * index + 2;

        if left < heap_len && array[left] > array[biggest] {
            biggest = left;
        }
        if right < heap_len && array[right] > array[biggest] {
            biggest = right;
        }

        if biggest != index {
            array.swap(index, biggest);
            sa_prefix.swap(index, biggest);
            index = biggest;
        } else {
            break;
        }
    }
}

fn build_heap(sa: &mut [usize], n_s: usize) {
    let mut i = n_s / 2;
    while i > 0 {
        let index = i - 1;
        heapify(sa, n_s, n_s, index);
        i -= 1;
    }
}

fn custom_binary_search(
    input: &[usize],
    array: &[usize],
    length: usize,
    value: usize,
    using_l_type: bool,
) -> usize {
    if length == 0 {
        return length;
    }
    let mut left = 0;
    let mut right = length - 1;

    let is_valid = |index: usize| {
        array[index] < length
            && (using_l_type || index == 0 || array[index.wrapping_sub(1)] != BH)
            && (!using_l_type || index == length - 1 || array[index + 1] != BH)
    };

    while right > left {
        let mid = left + (right - left) / 2;
        let mut checking = mid;
        while checking > left && !is_valid(checking) {
            checking -= 1;
        }

        let arr_val = if is_valid(checking) { array[checking] } else { usize::MAX };
        let inp_val = if arr_val < input.len() { input[arr_val] } else { usize::MAX };

        if is_valid(checking) && inp_val >= value {
            right = checking;
        } else {
            left = mid + 1;
        }
    }

    let arr_left = array[left];
    if left >= length || arr_left >= length || input[arr_left] != value {
        length
    } else {
        left
    }
}

fn custom_binary_search_last(
    input: &[usize],
    array: &[usize],
    length: usize,
    value: usize,
) -> usize {
    if length == 0 {
        return length;
    }
    let mut left = 0;
    let mut right = length - 1;

    let is_valid = |index: usize| {
        array[index] < length && (index == length - 1 || array[index + 1] != BH)
    };

    while right > left {
        let mid = left + (right - left) / 2 + ((right - left) % 2);
        let mut checking = mid;
        while checking < right && !is_valid(checking) {
            checking += 1;
        }

        let arr_val = if is_valid(checking) { array[checking] } else { usize::MAX };
        let inp_val = if arr_val < input.len() { input[arr_val] } else { usize::MAX };

        if is_valid(checking) && inp_val <= value {
            left = checking;
        } else {
            right = mid.saturating_sub(1);
        }
    }

    let arr_right = array[right];
    if right >= length || arr_right >= length || input[arr_right] != value {
        length
    } else {
        right
    }
}

// ----- Core Algorithm Functions -----

fn count_s_type(input: &[usize]) -> usize {
    let length = input.len();
    if length <= 1 {
        return 0;
    }
    let mut next_is_l = true;
    let mut count = 0;

    for i in 0..(length - 1) {
        let index = length - i - 2;
        let current_is_l =
            input[index] > input[index + 1] || (next_is_l && input[index] == input[index + 1]);
        if !current_is_l {
            count += 1;
        }
        next_is_l = current_is_l;
    }
    count
}

fn place_indices_of_type(input: &[usize], sa: &mut [usize], using_l_type: bool) {
    let length = sa.len();
    let mut next_is_l = false;
    let mut last_inserted_index = length;

    for i in 0..length {
        let index = length - i - 1;
        let current_is_l = index == length - 1
            || input[index] > input[index + 1]
            || (next_is_l && input[index] == input[index + 1]);

        if (current_is_l && using_l_type) || (!current_is_l && !using_l_type) {
            last_inserted_index = last_inserted_index.saturating_sub(1);
            if last_inserted_index < length {
                sa[last_inserted_index] = index;
            }
        }
        next_is_l = current_is_l;
    }
}

fn merge_sort_substrings(input: &[usize], sa: &mut [usize], n_s: usize, using_l_type: bool) {
    let length = sa.len();
    if n_s == 0 || length < n_s { return; }
    let mut step = 1;

    while step < length {
        let mut i = length - n_s;
        while i + step < length {
            let mid = i + step;
            let end = (i + step * 2).min(length);

            let (aux, target) = sa.split_at_mut(i);
            let left = &target[..mid - i];
            let right = &target[mid - i..end - i];

            let mut p1 = 0;
            let mut p2 = 0;
            let mut k = 0;

            while p1 < left.len() && p2 < right.len() {
                if compare_substrings(input, left[p1], right[p2], length, using_l_type)
                    == Ordering::Less
                {
                    aux[k] = left[p1];
                    p1 += 1;
                } else {
                    aux[k] = right[p2];
                    p2 += 1;
                }
                k += 1;
            }
            while p1 < left.len() {
                aux[k] = left[p1];
                p1 += 1;
                k += 1;
            }
            while p2 < right.len() {
                aux[k] = right[p2];
                p2 += 1;
                k += 1;
            }

            target[..end - i].copy_from_slice(&aux[..end - i]);
            i += step * 2;
        }
        step *= 2;
    }
}

fn construct_reduced_problem(
    input: &[usize],
    sa: &mut [usize],
    n_s: usize,
    using_l_type: bool,
) {
    let length = sa.len();
    if n_s == 0 || length < n_s { return; }
    let mut current_char = 0;
    sa[0] = current_char;

    let offset = length - n_s;
    for i in (offset + 1)..length {
        if compare_substrings(input, sa[i - 1], sa[i], length, using_l_type) != Ordering::Equal {
            current_char += 1;
        }
        sa[i - offset] = current_char;
    }
}

fn heap_sort_reduced_problem(sa: &mut [usize], n_s: usize) {
    if n_s == 0 { return; }
    build_heap(sa, n_s);
    let length = sa.len();

    let mut i = n_s - 1;
    while i > 0 {
        let (sa_prefix, sa_suffix) = sa.split_at_mut(length - n_s);
        sa_suffix.swap(0, i);
        sa_prefix.swap(0, i);

        heapify(sa, n_s, i, 0);
        i -= 1;
    }
}

fn restore_from_recursion(
    input: &[usize],
    sa: &mut [usize],
    n_s: usize,
    using_l_type: bool,
) {
    let length = sa.len();
    let mut next_is_l = false;
    let mut sum = 0;

    for i in 0..length {
        let index = length - i - 1;
        let current_is_l = index == length - 1
            || input[index] > input[index + 1]
            || (next_is_l && input[index] == input[index + 1]);

        if (!current_is_l && !using_l_type) || (current_is_l && using_l_type) {
            sum += 1;
            if n_s >= sum {
                sa[n_s - sum] = index;
            }
        }
        next_is_l = current_is_l;
    }

    let start = length.saturating_sub(n_s);
    for i in start..length {
        if sa[i] < sa.len() {
            sa[i] = sa[sa[i]];
        }
    }
}

fn preprocess(input: &[usize], sa: &mut [usize], n_s: usize, using_l_type: bool) {
    let length = sa.len();
    let mut next_is_l = false;
    let mut last_inserted_index = length.saturating_sub(n_s);

    for i in 0..length {
        let index = length - i - 1;
        let current_is_l = index == length - 1
            || input[index] > input[index + 1]
            || (next_is_l && input[index] == input[index + 1]);

        if (current_is_l && !using_l_type) || (!current_is_l && using_l_type) {
            if last_inserted_index > 0 {
                last_inserted_index -= 1;
                sa[last_inserted_index] = index;
            }
        }
        next_is_l = current_is_l;
    }

    // Must remain stable. We replaced unwrap_or with direct indexing for performance.
    sa.sort_by(|&a, &b| {
        let key_a = if a < input.len() { input[a] } else { usize::MAX };
        let key_b = if b < input.len() { input[b] } else { usize::MAX };
        key_a.cmp(&key_b)
    });
}

fn initialize_sa(input: &[usize], sa: &mut [usize], using_l_type: bool) {
    let length = sa.len();
    let mut next_is_l = false;

    for i in 0..length {
        let index = length - i - 1;
        let current_is_l = index == length - 1
            || input[index] > input[index + 1]
            || (next_is_l && input[index] == input[index + 1]);

        if (current_is_l && !using_l_type) || (!current_is_l && using_l_type) {
            let first = custom_binary_search(input, sa, length, input[index], using_l_type);
            let l = if using_l_type {
                custom_binary_search_last(input, sa, length, input[index])
            } else {
                first
            };

            if (using_l_type && l == 0) || (!using_l_type && l == length - 1) || l >= length || first >= length {
                next_is_l = current_is_l;
                continue;
            }

            let position = if using_l_type { l.wrapping_sub(1) } else { l + 1 };
            let further_position = if using_l_type { l.wrapping_sub(2) } else { l + 2 };

            if position < length && sa[position] >= length {
                next_is_l = current_is_l;
                continue;
            }

            let mut r_l = l;
            while r_l < length && sa[r_l] != index {
                if using_l_type {
                    if r_l == 0 {
                        r_l = usize::MAX;
                        break;
                    }
                    r_l -= 1;
                } else {
                    r_l += 1;
                }
            }

            if r_l >= length {
                next_is_l = current_is_l;
                continue;
            }

            let n_l = if using_l_type {
                if r_l >= first { r_l - first + 1 } else { 0 }
            } else {
                if r_l >= l { r_l - l + 1 } else { 0 }
            };

            if n_l == 0 {
                next_is_l = current_is_l;
                continue;
            }

            if using_l_type && first <= l && l < sa.len() {
                let slice_len = l + 1 - first;
                if slice_len > 0 {
                    let rot = n_l % slice_len;
                    sa[first..=l].rotate_left(rot);
                }
            }

            if n_l == 2 {
                if position < length {
                    sa[position] = BT;
                }
            } else if n_l == 3 {
                if position < length {
                    sa[position] = BH;
                }
                if further_position < length {
                    sa[further_position] = BT;
                }
            } else if n_l > 3 {
                if position < length {
                    sa[position] = BH;
                }
                if further_position < length {
                    sa[further_position] = E;
                }
                if using_l_type {
                    if l + 1 >= n_l {
                        let idx = l + 1 - n_l;
                        if idx < length {
                            sa[idx] = BT;
                        }
                    }
                } else {
                    let idx = l + n_l - 1;
                    if idx < length {
                        sa[idx] = BT;
                    }
                }
            }
        }
        next_is_l = current_is_l;
    }
}

// ----- Context for Induced Sorting Caches -----

struct Case4Cache {
    l_cached: usize,
    s_cached: usize,
}

impl Case4Cache {
    #[inline]
    fn new() -> Self {
        Self {
            l_cached: usize::MAX,
            s_cached: usize::MAX,
        }
    }

    fn case4(&mut self, input: &[usize], sa: &mut [usize], l: usize, j: usize) {
        let length = sa.len();
        if self.l_cached != usize::MAX && l == self.l_cached {
            return;
        } else {
            self.l_cached = l;
        }

        let mut r_l = l + 1;
        while r_l < length
            && sa[r_l] != R1
            && sa[r_l] < length
            && input[sa[r_l]] == input[j]
        {
            r_l += 1;
        }
        if r_l < length
            && (sa[r_l] == R1 || (sa[r_l] < length && input[sa[r_l]] == input[j]))
        {
            sa[r_l] = j;
        }
    }

    fn case4_s(&mut self, input: &[usize], sa: &mut [usize], l: usize, j: usize) {
        if self.s_cached != usize::MAX && l == self.s_cached {
            return;
        } else {
            self.s_cached = l;
        }
        let length = sa.len();

        let mut r_l = l;
        while r_l > 0
            && sa[r_l - 1] != R1
            && sa[r_l - 1] < length
            && input[sa[r_l - 1]] == input[j]
        {
            r_l -= 1;
        }
        if r_l != 0 && sa[r_l - 1] == R1 {
            sa[r_l - 1] = j;
        }
    }

    fn sort_s_body(&mut self, input: &[usize], sa: &mut [usize], j: usize) {
        let length = sa.len();
        if j == length - 1 || input[j] > input[j + 1] {
            return;
        }

        let l = custom_binary_search_last(input, sa, length, input[j]);
        if l == 0 || l >= length {
            return;
        } else if input[j] == input[j + 1] {
            if !(sa[l - 1] == BH || sa[l - 1] == R2 || sa[l - 1] == R1 || sa[l - 1] == BT) {
                self.case4_s(input, sa, l, j);
                return;
            }
        }

        let prev = sa[l - 1];

        if prev == BT {
            sa[l] = j;
            sa[l - 1] = R1;
        } else if prev == R1 {
            sa[l - 1] = j;
        } else if prev == BH {
            if l < 2 {
                return;
            }
            let l_minus_2 = l - 2;
            if sa[l_minus_2] == E {
                sa[l] = j;
                sa[l_minus_2] = 1;
            } else if sa[l_minus_2] == BT {
                sa[l] = j;
                sa[l_minus_2] = R2;
            } else if sa[l_minus_2] == R2 {
                sa[l - 1] = j;
                sa[l_minus_2] = R1;
            } else if sa[l_minus_2] == R1 {
                sa[l_minus_2] = j;
            } else {
                let c = sa[l_minus_2];
                if l < c + 2 {
                    return;
                }
                let target_idx = l - c - 2;
                if sa[target_idx] != BT {
                    sa[target_idx] = j;
                    sa[l_minus_2] = c + 1;
                } else {
                    let r_l = l - 2 - c;
                    let src_start = r_l + 1;
                    let src_end = l - 2;
                    let dest = r_l + 2;
                    if src_start <= src_end && dest + (src_end - src_start) <= length {
                        sa.copy_within(src_start..src_end, dest);
                        sa[r_l + 1] = j;
                        sa[l - 1] = R2;
                    }
                }
            }
        } else if prev == R2 {
            if l < 2 {
                return;
            }
            let mut r_l = l - 2;
            while r_l < length && sa[r_l] != BT {
                if r_l == 0 {
                    return;
                }
                r_l -= 1;
            }
            let src_start = r_l + 1;
            let src_end = l - 1;
            let dest = r_l + 2;
            if src_start <= src_end && dest + (src_end - src_start) <= length {
                sa.copy_within(src_start..src_end, dest);
                sa[r_l + 1] = j;
                sa[r_l] = R1;
            }
        } else {
            self.case4_s(input, sa, l, j);
        }
    }

    fn sort_l_body(&mut self, input: &[usize], sa: &mut [usize], j: usize) {
        let length = sa.len();
        if j != length - 1 && input[j] < input[j + 1] {
            return;
        }

        let l = custom_binary_search(input, sa, length, input[j], false);
        if l >= length - 1 {
            return;
        } else if j != length - 1 && input[j] == input[j + 1] {
            if !(sa[l + 1] == BH || sa[l + 1] == R2 || sa[l + 1] == R1 || sa[l + 1] == BT) {
                self.case4(input, sa, l, j);
                return;
            }
        }

        if sa[l + 1] == BT {
            sa[l] = j;
            sa[l + 1] = R1;
        } else if sa[l + 1] == R1 {
            sa[l + 1] = j;
        } else if sa[l + 1] == BH {
            if l + 2 >= length {
                return;
            }
            if sa[l + 2] == E {
                sa[l] = j;
                sa[l + 2] = 1;
            } else if sa[l + 2] == BT {
                sa[l] = j;
                sa[l + 2] = R2;
            } else if sa[l + 2] == R2 {
                sa[l + 1] = j;
                sa[l + 2] = R1;
            } else if sa[l + 2] == R1 {
                sa[l + 2] = j;
            } else {
                let c = sa[l + 2];
                if l + c + 2 < length && sa[l + c + 2] != BT {
                    sa[l + c + 2] = j;
                    sa[l + 2] = c + 1;
                } else {
                    let r_l = l + 2 + c;
                    let src_start = l + 3;
                    let src_end = r_l;
                    let dest = l + 2;
                    if src_start <= src_end && src_end <= length && dest + (src_end - src_start) <= length {
                        sa.copy_within(src_start..src_end, dest);
                        if r_l > 0 {
                            sa[r_l - 1] = j;
                        }
                        sa[l + 1] = R2;
                    }
                }
            }
        } else if sa[l + 1] == R2 {
            let mut r_l = l + 2;
            while r_l < length - 1 && sa[r_l] != BT {
                r_l += 1;
            }
            let src_start = l + 2;
            let src_end = r_l;
            let dest = l + 1;
            if src_start <= src_end && src_end <= length && dest + (src_end - src_start) <= length {
                sa.copy_within(src_start..src_end, dest);
                if r_l > 0 {
                    sa[r_l - 1] = j;
                }
                sa[r_l] = R1;
            }
        } else {
            self.case4(input, sa, l, j);
        }
    }
}

fn sort_s(input: &[usize], sa: &mut [usize], cache: &mut Case4Cache) {
    let length = sa.len();
    let mut iter = 1;
    while iter <= length {
        let i = length - iter;
        if sa[i] == 0 {
            iter += 1;
            continue;
        }

        if sa[i] >= length {
            if sa[i] == BH {
                iter += 1;
            }
            iter += 1;
            continue;
        }

        let mut j = sa[i].wrapping_sub(1);
        if j < length {
            cache.sort_s_body(input, sa, j);

            while sa[i] != 0 && sa[i] < length && sa[i].wrapping_sub(1) != j {
                j = sa[i].wrapping_sub(1);
                if j >= length {
                    break;
                }
                cache.sort_s_body(input, sa, j);
            }
        }

        iter += 1;
    }
}

fn sort_l(input: &[usize], sa: &mut [usize], cache: &mut Case4Cache) {
    let length = sa.len();
    if length == 0 { return; }
    cache.sort_l_body(input, sa, length - 1);

    let mut i = 0;
    while i < length {
        if sa[i] == 0 {
            i += 1;
            continue;
        }

        if sa[i] >= length {
            if sa[i] == BH {
                i += 1;
            }
            i += 1;
            continue;
        }

        let mut j = sa[i].wrapping_sub(1);
        if j < length {
            cache.sort_l_body(input, sa, j);

            while sa[i] != 0 && sa[i] < length && sa[i].wrapping_sub(1) != j {
                j = sa[i].wrapping_sub(1);
                if j >= length {
                    break;
                }
                cache.sort_l_body(input, sa, j);
            }
        }

        i += 1;
    }
}

fn induced_sorting(input: &[usize], sa: &mut [usize], using_l_type: bool) {
    let mut cache = Case4Cache::new();
    if using_l_type {
        sort_s(input, sa, &mut cache);
    } else {
        sort_l(input, sa, &mut cache);
    }
}

// ----- Main Entry Points -----

pub fn optimal_suffix_array_first(
    input: &[usize],
    sa: &mut [usize],
    prev_n_s: &mut usize,
) -> bool {
    let length = sa.len();
    if length < 2 {
        if length == 1 {
            sa[0] = 0;
        }
        return false;
    }

    let mut n_s = count_s_type(input);
    if n_s == 0 {
        for i in 0..length {
            sa[i] = length - i - 1;
        }
        return false;
    }

    let using_l_type = (length - n_s) < n_s;
    if using_l_type {
        n_s = length - n_s;
    }
    *prev_n_s = n_s;

    // Step 1
    place_indices_of_type(input, sa, using_l_type);

    // Step 2
    merge_sort_substrings(input, sa, n_s, using_l_type);

    // Step 3
    construct_reduced_problem(input, sa, n_s, using_l_type);

    // Step 4
    heap_sort_reduced_problem(sa, n_s);

    true
}

pub fn optimal_suffix_array_second(input: &[usize], sa: &mut [usize]) {
    let length = sa.len();
    let mut n_s = count_s_type(input);

    let using_l_type = (length - n_s) < n_s;
    if using_l_type {
        n_s = length - n_s;
    }

    // Step 5
    restore_from_recursion(input, sa, n_s, using_l_type);

    // Step 6
    preprocess(input, sa, n_s, using_l_type);

    // Step 7
    initialize_sa(input, sa, using_l_type);

    // Step 8
    induced_sorting(input, sa, using_l_type);
}

unsafe fn optimal_suffix_array(
    input: *mut usize,
    sa: *mut usize,
    length: usize,
) {
    let mut prev_n_s: usize = 0;
    
    let input_slice = std::slice::from_raw_parts(input, length);
    let sa_slice = std::slice::from_raw_parts_mut(sa, length);
    
    let mut needs_recursion = optimal_suffix_array_first(input_slice, sa_slice, &mut prev_n_s);
    if !needs_recursion {
        return;
    }

    let mut old_output = sa;
    let mut old_length = length;
    let mut new_output_index = 0;
    let mut old_old_length = old_length;
    let mut recursion_depth = 0;
    let mut new_input = std::ptr::null_mut();

    while needs_recursion {
        new_input = old_output;
        let new_length = prev_n_s;
        let new_output = old_output.add(old_length - new_length);
        
        let prev_had_space = (*new_input) & MSB != 0;
        if prev_had_space {
            *new_input &= !MSB;
        }

        let new_input_slice = std::slice::from_raw_parts(new_input, new_length);
        let new_output_slice = std::slice::from_raw_parts_mut(new_output, new_length);
        needs_recursion = optimal_suffix_array_first(new_input_slice, new_output_slice, &mut prev_n_s);
        
        if prev_had_space {
            *new_input |= MSB;
        }

        if needs_recursion {
            new_output_index += old_length - new_length;
            recursion_depth += 1;

            if new_length * 2 != old_length {
                *new_output |= MSB;
                *new_output.sub(1) = old_length;
            }

            old_output = new_output;
        }

        old_old_length = old_length;
        old_length = new_length;
    }
    
    let mut old_input = new_input;

    while recursion_depth > 0 {
        recursion_depth -= 1;

        old_length = old_old_length;
        let had_space = (*old_input) & MSB != 0;
        if had_space {
            *old_input &= !MSB;
            old_old_length = *old_input.sub(1);
        } else {
            old_old_length *= 2;
        }
        old_output = old_input;

        new_output_index -= old_old_length - old_length;
        old_input = sa.add(new_output_index);

        let current_have_space = (*old_input) & MSB != 0;
        if current_have_space {
            *old_input &= !MSB;
        }

        let old_input_slice = std::slice::from_raw_parts(old_input, old_length);
        let old_output_slice = std::slice::from_raw_parts_mut(old_output, old_length);
        optimal_suffix_array_second(old_input_slice, old_output_slice);

        if current_have_space {
            *old_input |= MSB;
        }
    }
    
    let input_slice = std::slice::from_raw_parts(input, length);
    let sa_slice = std::slice::from_raw_parts_mut(sa, length);
    optimal_suffix_array_second(input_slice, sa_slice);
}

/// Convenience wrapper for byte slice or string slice inputs
pub fn optimal_suffix_array_bytes(input_chars: &[u8], sa: &mut [usize]) {
    let mut input = chars_to_usize_array(input_chars);
    unsafe { 
        optimal_suffix_array(input.as_mut_ptr(), sa.as_mut_ptr(), input.len()) 
    };
}