use crate::saca::optimal_suffix_array_bytes;

pub fn bwt(s: &[u8]) -> Vec<u8> {
    let len = s.len();
    let mut sa = vec![0; len];
    optimal_suffix_array_bytes(&s, &mut sa);

    let mut bbwt = Vec::with_capacity(len); // todo can re-use sa
    for &i in &sa {
        let prev_index = if i == 0 { s.len() - 1 } else { i - 1 };
        bbwt.push(s[prev_index]);
    }
    bbwt
}