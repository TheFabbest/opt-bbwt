// todo since that only happens for periods, can't I just do start = start % d and find the first msp?

pub fn from_period_slice(s: &[u8], start: usize, period_length: usize) -> Vec<u8> {
    let start = start % period_length;
    let word: Vec<u8> = if start + period_length <= s.len() {
        s[start..start + period_length].to_vec()
    } else {
        // this should never happen, warn.
        eprintln!("This should never happen!");
        [&s[start..], &s[..start + period_length - s.len()]].concat()
    };
    word
}

pub fn rotate_string(s: &[u8], start: usize) -> Vec<u8> {
    let len = s.len();
    let start = start % len;
    let rotated: Vec<u8> = [&s[start..], &s[..start]].concat();
    rotated
}