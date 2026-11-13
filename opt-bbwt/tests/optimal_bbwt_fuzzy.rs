use opt_bbwt::{
    optimal_bbwt::{chain::optimal_bbwt_chain, naive::naive_optimal_bbwt},
    transform::bbwt::bbwt,
};

fn rotate(text: &[u8], shift: usize) -> Vec<u8> {
    [text[shift..].as_ref(), &text[..shift]].concat()
}

fn bbwt_runs(text: &[u8]) -> usize {
    if text.is_empty() {
        return 0;
    }
    bbwt(text)
        .windows(2)
        .filter(|window| window[0] != window[1])
        .count()
        + 1
}

fn assert_chain_is_optimal(input: &[u8]) {
    let chain_shift = optimal_bbwt_chain(input, true);
    let naive_shift = naive_optimal_bbwt(input);

    assert_eq!(
        bbwt_runs(&rotate(input, chain_shift)),
        bbwt_runs(&rotate(input, naive_shift)),
        "chain result is not optimal for input {input:?}"
    );
}

#[test]
fn chain_matches_naive_on_fuzzy_inputs() {
    let mut state = 0x4d595df4d0f33173u64;

    for case in 0..512 {
        let length = (next_random(&mut state) % 65) as usize;
        let mut input = vec![0; length];
        for byte in &mut input {
            *byte = b'a' + (next_random(&mut state) % 4) as u8;
        }

        assert_chain_is_optimal(&input);
        assert!(
            optimal_bbwt_chain(&input, false) < input.len().max(1),
            "chain returned an invalid shift in case {case} for input {input:?}"
        );
    }
}

#[test]
fn chain_matches_naive_for_all_binary_strings_up_to_length_ten() {
    for length in 0..=10 {
        for mask in 0..(1usize << length) {
            let input: Vec<u8> = (0..length)
                .map(|index| if mask & (1 << index) == 0 { b'a' } else { b'b' })
                .collect();
            assert_chain_is_optimal(&input);
        }
    }
}

fn next_random(state: &mut u64) -> u64 {
    *state ^= *state << 7;
    *state ^= *state >> 9;
    *state
}
