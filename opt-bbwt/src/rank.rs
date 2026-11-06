/// Converts suffix-array positions into ranks.
pub fn sa_to_rank(sa: &[usize]) -> Vec<usize> {
    let n = sa.len();
    let mut rank = vec![0; n];
    for (i, &suffix_index) in sa.iter().enumerate() {
        rank[suffix_index] = i;
    }
    rank
}
