use opt_bbwt::lyndon::{self, CanonicalRotationFinder, conjugatesfactors::for_each_factor_of_each_conjugate, duval, shiloach};

fn main() {
    let s = b"bbaabca";
    //lyndon::print_all_factors(s);

    // text is a randomly generated DNA string (&[u8]) of length 1000
    let generated: Vec<u8> = (0..100).map(|_| {
        let n = rand::random::<u8>() % 4;
        match n {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }).collect();
    let text = generated.as_slice();

    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(text);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| text[(start + k) % text.len()]).collect();


    // bench and compare
    let start_time = std::time::Instant::now();
    let naive_result = opt_bbwt::optimalbbwt::naive_optimal_bbwt(aperiodic_least_rotation.as_slice());
    let naive_time = start_time.elapsed();
    println!("Naive optimal BBWT index: {}. Time taken: {:?}", naive_result, naive_time);

    let start_time = std::time::Instant::now();
    let optimal_result = opt_bbwt::optimalbbwt::sort_all_factors(aperiodic_least_rotation.as_slice());
    let optimal_time = start_time.elapsed();
    println!("Optimal BBWT index: {}. Time taken: {:?}", optimal_result, optimal_time);
}