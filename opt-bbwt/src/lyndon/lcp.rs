/// Computes the longest-common-prefix array from `text` and its suffix array.
pub fn lcp(text: &[u8], sa: &[usize], rank: &[usize]) -> Vec<usize> {
    let n = text.len();
    let mut lcp = vec![0; n];
    let mut h = 0;

    for i in 0..n {
        if rank[i] > 0 {
            let j = sa[rank[i] - 1];
            while i + h < n && j + h < n && text[i + h] == text[j + h] {
                h += 1;
            }
            lcp[rank[i]] = h;
            h = h.saturating_sub(1);
        }
    }

    lcp
}
