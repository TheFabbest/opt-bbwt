use crate::lyndon::shiloach;
use crate::transform::bbwt::bbwt;

pub fn naive_optimal_bbwt(text: &[u8]) -> usize {
    if text.is_empty() {
        return 0;
    }
    let (start, n) = shiloach(text);

    let s = (0..n)
        .map(|k| text[(start + k) % text.len()])
        .collect::<Vec<u8>>();

    let mut min_runs_index = 0;
    let mut min_runs = usize::MAX;
    for i in 0..n {
        let conjugate: Vec<u8> = s[i..].iter().chain(&s[..i]).copied().collect();
        let bbwt_result = bbwt(&conjugate);
        let runs = bbwt_result.windows(2).filter(|w| w[0] != w[1]).count() + 1; // count runs in BBWT
        if runs < min_runs {
            min_runs = runs;
            min_runs_index = i;
        }
    }

    min_runs_index
}
