use opt_bbwt::{
    optimalbbwt,
    suffix_array::circular_suffix_array,
    transform,
};
use opt_bbwt::optimalbbwt::{
    globalsort::optimal_bbwt_global_sort,
    naive::naive_optimal_bbwt,
    optimalbbwt::online_bbwt_check,
    reuse::optimal_bbwt_reuse,
};
use opt_bbwt::suffix_array::{
    circular_sais::circular_suffix_array as circular_sais,
    naive_circular::naive_circular_suffix_array,
    optimal::optimal_suffix_array_bytes,
    sais::sais,
};
use opt_bbwt::transform::bbwt::naive::naive_bbwt;

#[test]
fn suffix_array_implementations_agree() {
    for input in [b"abracadabra".as_slice(), b"banana", b"aaaa"] {
        let optimal = optimal_suffix_array_bytes(input);
        let sais_output = sais(input);
        assert_eq!(optimal, sais_output);
        assert_eq!(sais(input), sais_output);
    }
}

#[test]
fn public_default_and_variant_entry_points_are_consistent() {
    let input = b"banana";
    assert_eq!(transform::bbwt::bbwt(input), naive_bbwt(input));
    assert_eq!(transform::bwt::bwt(input), {
        let mut output = vec![0; input.len()];
        transform::bwt::bwt_into(input, &mut output);
        output
    });
    assert_eq!(
        optimalbbwt::optimal_bbwt(input),
        online_bbwt_check(input)
    );
}

#[test]
fn optimal_bbwt_variants_accept_empty_input() {
    assert_eq!(optimalbbwt::optimal_bbwt(b""), 0);
    assert_eq!(naive_optimal_bbwt(b""), 0);
    assert_eq!(optimal_bbwt_global_sort(b""), 0);
    assert_eq!(optimal_bbwt_reuse(b""), 0);
    assert_eq!(online_bbwt_check(b""), 0);
}

#[test]
fn transform_implementations_agree() {
    for input in [b"banana".as_slice(), b"abracadabra", b""] {
        assert_eq!(transform::bbwt::bbwt(input), naive_bbwt(input));

        let mut bbwt_output = vec![0; input.len()];
        transform::bbwt::bbwt_into(input, &mut bbwt_output);
        assert_eq!(transform::bbwt::bbwt(input), bbwt_output);

        let mut bwt_output = vec![0; input.len()];
        transform::bwt::bwt_into(input, &mut bwt_output);
        assert_eq!(transform::bwt::bwt(input), bwt_output);
    }
}

#[test]
fn circular_suffix_array_implementations_agree() {
    for input in [b"banana".as_slice(), b"abracadabra", b"aaaa"] {
        let output = circular_suffix_array(input);
        assert_eq!(circular_sais(input), output);
        assert_eq!(output, naive_circular_suffix_array(input));
    }
}

#[test]
fn inverse_bbwt_into_matches_allocating_api() {
    let input = b"banana";
    let transformed = transform::bbwt::bbwt(input);
    let mut output = vec![0; input.len()];
    transform::bbwt::unbbwt_into(&transformed, &mut output);
    assert_eq!(transform::bbwt::unbbwt(&transformed), output);
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
        let result_through_w = optimal_bbwt_global_sort(input);
        let result_reuse = optimal_bbwt_reuse(input);
        let result_online = online_bbwt_check(input);
        let naive_result = naive_optimal_bbwt(input);

        let rotation_through_w = [
            input[result_through_w..].as_ref(),
            &input[..result_through_w],
        ]
        .concat();
        let rotation_reuse = [input[result_reuse..].as_ref(), &input[..result_reuse]].concat();
        let rotation_online = [input[result_online..].as_ref(), &input[..result_online]].concat();
        let naive_rotation = [input[naive_result..].as_ref(), &input[..naive_result]].concat();
        let runs = |text: &[u8]| {
            transform::bbwt::bbwt(text)
                .windows(2)
                .filter(|w| w[0] != w[1])
                .count()
                + 1
        };
        assert_eq!(runs(&rotation_through_w), runs(&naive_rotation));
        assert_eq!(runs(&rotation_reuse), runs(&naive_rotation));
        assert_eq!(runs(&rotation_online), runs(&naive_rotation));
    }
}
