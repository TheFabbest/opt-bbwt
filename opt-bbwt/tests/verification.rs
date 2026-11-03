use opt_bbwt::{
    bbwt::{bbwt, bbwt_into, unbbwt, unbbwt_into},
    bwt::{bwt, bwt_into},
    naive,
    optimalbbwt::{naive_optimal_bbwt, optimal_bbwt_through_w},
    saca::{
        circular_suffix_array, circular_suffix_array_into, optimal_suffix_array_bytes, sais,
        sais_into,
    },
};

#[test]
fn suffix_array_implementations_agree() {
    for input in [b"abracadabra".as_slice(), b"banana", b"aaaa"] {
        let mut optimal = vec![0; input.len()];
        let mut sais_output = vec![0; input.len()];
        optimal_suffix_array_bytes(input, &mut optimal);
        sais_into(input, &mut sais_output);
        assert_eq!(optimal, sais_output);
        assert_eq!(sais(input), sais_output);
    }
}

#[test]
fn transform_implementations_agree() {
    for input in [b"banana".as_slice(), b"abracadabra", b""] {
        assert_eq!(bbwt(input), naive::bbwt::bbwt(input));

        let mut bbwt_output = vec![0; input.len()];
        bbwt_into(input, &mut bbwt_output);
        assert_eq!(bbwt(input), bbwt_output);

        let mut bwt_output = vec![0; input.len()];
        bwt_into(input, &mut bwt_output);
        assert_eq!(bwt(input), bwt_output);
    }
}

#[test]
fn circular_suffix_array_implementations_agree() {
    for input in [b"banana".as_slice(), b"abracadabra", b"aaaa"] {
        let mut output = vec![0; input.len()];
        circular_suffix_array_into(input, &mut output);
        assert_eq!(circular_suffix_array(input), output);
        assert_eq!(output, naive::csa::circular_suffix_array_naive(input));
    }
}

#[test]
fn inverse_bbwt_into_matches_allocating_api() {
    let input = b"banana";
    let transformed = bbwt(input);
    let mut output = vec![0; input.len()];
    unbbwt_into(&transformed, &mut output);
    assert_eq!(unbbwt(&transformed), output);
    assert_eq!(output, input);
}

#[test]
fn optimal_bbwt_matches_quadratic_reference() {
    for input in [
        b"abcabcabc".as_slice(),
        b"aaaaaa",
        b"banana",
        b"mississippi",
    ] {
        let result = optimal_bbwt_through_w(input);
        let naive_result = naive_optimal_bbwt(input);
        let rotation = [input[result..].as_ref(), &input[..result]].concat();
        let naive_rotation = [input[naive_result..].as_ref(), &input[..naive_result]].concat();
        let runs = |text: &[u8]| bbwt(text).windows(2).filter(|w| w[0] != w[1]).count() + 1;
        assert_eq!(runs(&rotation), runs(&naive_rotation));
    }
}
