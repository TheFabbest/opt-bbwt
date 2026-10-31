use opt_bbwt::lyndon::{CanonicalRotationFinder, shiloach};

fn main() {
    //let s = b"bbaabca";
    //lyndon::print_all_factors(s);

    // text is the first characters from the dna corpus
    let corpus_size = 3000;
    let corpus = opt_bbwt::file::load_corpus("dna", corpus_size);
    let text = corpus.expect("Failed to load corpus");

    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(&text);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| text[(start + k) % text.len()]).collect();


    // bench and compare
    let start_time = std::time::Instant::now();
    let naive_result = opt_bbwt::optimalbbwt::naive_optimal_bbwt(aperiodic_least_rotation.as_slice());
    let naive_time = start_time.elapsed();
    println!("Naive optimal BBWT index: {}. Time taken: {:?}", naive_result, naive_time);

    // let start_time = std::time::Instant::now();
    // let mine = opt_bbwt::optimalbbwt::sort_all_factors(aperiodic_least_rotation.as_slice());
    // let mine_time = start_time.elapsed();\
    // println!("My optimal BBWT index: {}. Time taken: {:?}", mine, mine_time);

    let start_time = std::time::Instant::now();
    let runs = opt_bbwt::optimalbbwt::optimal_bbwt_through_w(&aperiodic_least_rotation);
    let optimal_time = start_time.elapsed();
    println!("Optimal BBWT index (through w): {}. Time taken: {:?}", runs, optimal_time);
}