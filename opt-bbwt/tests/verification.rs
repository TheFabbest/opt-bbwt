use opt_bbwt::{
    bbwt::bbwt,
    naive,
    saca::{optimal_suffix_array_bytes, sais},
};

#[test]
fn optimized_suffix_array_matches_sais() {
    let input = b"abracadabra";
    let mut optimized = vec![0; input.len()];
    optimal_suffix_array_bytes(input, &mut optimized);
    assert_eq!(optimized, sais(input));
}

#[test]
fn optimized_bbwt_matches_reference() {
    for input in [b"banana".as_slice(), b"abracadabra", b""] {
        assert_eq!(bbwt(input), naive::bbwt::bbwt(input));
    }
}
