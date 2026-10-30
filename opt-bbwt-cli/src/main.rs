use opt_bbwt::lyndon::{self, CanonicalRotationFinder, conjugatesfactors::for_each_factor_of_each_conjugate, duval, shiloach};

fn main() {
    let text = b"bbaabca";
    //lyndon::print_all_factors(text);
    
    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(text);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| text[(start + k) % text.len()]).collect();

    opt_bbwt::optimalbbwt::sort_all_factors(aperiodic_least_rotation.as_slice());

}