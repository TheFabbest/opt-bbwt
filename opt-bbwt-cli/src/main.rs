use opt_bbwt::{bbwt::bbwt, lyndon::{CanonicalRotationFinder, lyndon_suffix_table, lyndon_table, shiloach}, rank::sa_to_rank, saca::{optimal_suffix_array_bytes}, string::rotate_string};
use opt_bbwt::heuristics::{most_factors, most_repeated_factors, most_repeated_factors_propotional};

fn main() {
    // text is the first characters from the dna corpus
    let corpus_size = 5000000;
    let corpus = opt_bbwt::file::load_corpus("dna", corpus_size);
    let text = corpus.expect("Failed to load corpus");

    let finder = shiloach::Shiloach;
    let (start, period_length) = finder.find(&text);

    // first `period_length` characters of the least rotation, wrapping around the end
    let aperiodic_least_rotation: Vec<u8> =
        (0..period_length).map(|k| text[(start + k) % text.len()]).collect();


    // bench and compare
    // let start_time = std::time::Instant::now();
    // let best_index = opt_bbwt::optimalbbwt::optimal_bbwt_through_w(&aperiodic_least_rotation);
    // let optimal_time = start_time.elapsed();
    // println!("Optimal BBWT index (through w): {}. Time taken: {:?}", best_index, optimal_time);
    // let through_w_best_runs = bbwt(&rotate_string(&aperiodic_least_rotation, best_index)).windows(2).filter(|w| w[0] != w[1]).count();

    // // let start_time = std::time::Instant::now();
    // // let naive_result = opt_bbwt::optimalbbwt::naive_optimal_bbwt(aperiodic_least_rotation.as_slice());
    // // let naive_time = start_time.elapsed();
    // // println!("Naive optimal BBWT index: {}. Time taken: {:?}", naive_result, naive_time);
    // // let naive_best_runs = bbwt(&rotate_string(&aperiodic_least_rotation, naive_result)).windows(2).filter(|w| w[0] != w[1]).count();

    // // assert_eq!(naive_best_runs, through_w_best_runs, "Naive best runs: {}, Through W best runs: {}", naive_best_runs, through_w_best_runs);

    // // let mut worst_runs = 0;
    // // for i in 0..aperiodic_least_rotation.len() {
    // //     let rotated = rotate_string(&aperiodic_least_rotation, i);
    // //     let runs = bbwt(&rotated).windows(2).filter(|w| w[0] != w[1]).count();
    // //     if runs > worst_runs {
    // //         worst_runs = runs;
    // //     }
    // // }

    // //println!("Worst case runs: {}", worst_runs);
    // println!("Best number of runs: {}", through_w_best_runs);

    // use heuristics
    let start_time = std::time::Instant::now();
    let mut sa = vec![0; aperiodic_least_rotation.len()];
    optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);
    let rank = sa_to_rank(&sa);
    let lyn = &lyndon_table(&aperiodic_least_rotation, &rank);
    let lyn_s = &lyndon_suffix_table(&aperiodic_least_rotation);
    let mrp = most_repeated_factors(&aperiodic_least_rotation, lyn, lyn_s);
    let heuristic_time = start_time.elapsed();
    let runs = bbwt(&rotate_string(&aperiodic_least_rotation, mrp)).windows(2).filter(|w| w[0] != w[1]).count();
    println!("heuristic suggested: {}, that has {} runs. Time taken: {:?}", mrp, runs, heuristic_time);

    let start_time = std::time::Instant::now();
    let mut sa = vec![0; aperiodic_least_rotation.len()];
    optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);
    let rank: Vec<usize> = sa_to_rank(&sa);
    let lyn = &lyndon_table(&aperiodic_least_rotation, &rank);
    let lyn_s = &lyndon_suffix_table(&aperiodic_least_rotation);
    let mrp = most_factors(lyn, lyn_s);
    let heuristic_time = start_time.elapsed();
    let runs = bbwt(&rotate_string(&aperiodic_least_rotation, mrp)).windows(2).filter(|w| w[0] != w[1]).count();
    println!("heuristic suggested: {}, that has {} runs. Time taken: {:?}", mrp, runs, heuristic_time);

    let start_time = std::time::Instant::now();
    let mut sa = vec![0; aperiodic_least_rotation.len()];
    optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut sa);
    let rank = sa_to_rank(&sa);
    let lyn = &lyndon_table(&aperiodic_least_rotation, &rank);
    let lyn_s = &lyndon_suffix_table(&aperiodic_least_rotation);
    let mrp = most_repeated_factors_propotional(&aperiodic_least_rotation, lyn, lyn_s);
    let heuristic_time = start_time.elapsed();
    let runs = bbwt(&rotate_string(&aperiodic_least_rotation, mrp)).windows(2).filter(|w| w[0] != w[1]).count();
    println!("heuristic suggested: {}, that has {} runs. Time taken: {:?}", mrp, runs, heuristic_time);


    // normal bbwt
    if text.len() != aperiodic_least_rotation.len() {
        println!("Warning: The text length is not equal to the aperiodic least rotation length. This may affect the following runs count.");
    }
    let start = std::time::Instant::now();
    let bbwt = bbwt(&text);
    let bbwt_time = start.elapsed();
    println!("Normal BBWT time taken: {:?}", bbwt_time);
    let runs = bbwt.windows(2).filter(|w| w[0] != w[1]).count();
    println!("Normal BBWT runs: {}", runs);

    // normal bwt
    let start = std::time::Instant::now();
    let bwt = opt_bbwt::bwt::bwt(&text);
    let bwt_time = start.elapsed();
    println!("Normal BWT time taken: {:?}", bwt_time);
    let runs = bwt.windows(2).filter(|w| w[0] != w[1]).count();
    println!("Normal BWT runs: {}", runs);
}