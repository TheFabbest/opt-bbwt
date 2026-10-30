/// Computes the starting indices of the Lyndon factors of `s` using Duval's algorithm in O(n) time and O(1) extra space.
pub fn duval(s: &[u8]) -> Vec<usize> {
    let n = s.len();
    let mut result = Vec::new();
    let mut i = 0;

    while i < n {
        let mut j = i + 1;
        let mut k = i;

        while j < n && s[k] <= s[j] {
            if s[k] < s[j] {
                k = i;
            } else {
                k += 1;
            }
            j += 1;
        }

        while i <= k {
            result.push(i);
            i += j - k;
        }
    }

    result
}
