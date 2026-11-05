use clap::Parser;
use opt_bbwt::{
    bbwt::bbwt, heuristics::{most_factors, most_repeated_factors, most_repeated_factors_proportional}, lyndon::{lyndon_suffix_table, lyndon_table, shiloach}, optimalbbwt::optimalbbwt, rank::sa_to_rank, saca::{optimal_suffix_array_bytes, sais}, string::rotate_string,
};

#[derive(Debug, Parser)]
#[command(about = "Run BBWT and heuristic measurements on a corpus")]
struct Args {
    /// Name of the file under ./corpora.
    #[arg(short, long, default_value = "dna")]
    corpus: String,

    /// Maximum number of bytes to read from the corpus.
    #[arg(short, long, default_value_t = 50_000)]
    length: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let text = opt_bbwt::file::load_corpus(&args.corpus, args.length)?;
    if text.is_empty() {
        return Err("the selected corpus is empty".into());
    }

    let (start, period_length) = shiloach(&text);
    let aperiodic_least_rotation: Vec<u8> = (0..period_length)
        .map(|k| text[(start + k) % text.len()])
        .collect();

    // benchmark online
    let start = std::time::Instant::now();
    let rotation_online = optimalbbwt::optimal_bbwt_chain(&aperiodic_least_rotation, true);
    let online_duration = start.elapsed();
    println!("optimal_bbwt_online: rotation={rotation_online}, time={online_duration:?}");
    return Ok(());
    
    // compare optimal bbwt implementations
    let start = std::time::Instant::now();
    let naive_result = optimalbbwt::naive_optimal_bbwt(&aperiodic_least_rotation);
    let naive_duration = start.elapsed();
    let rotation_through_w = optimalbbwt::optimal_bbwt_through_w(&aperiodic_least_rotation);
    let through_w_duration = start.elapsed() - naive_duration;
    let rotation_reuse = optimalbbwt::optimal_bbwt_reuse(&aperiodic_least_rotation);
    let reuse_duration = start.elapsed() - naive_duration - through_w_duration;
    let rotation_online = optimalbbwt::optimal_bbwt_chain(&aperiodic_least_rotation, true);
    let online_duration = start.elapsed() - naive_duration - through_w_duration - reuse_duration;

    println!("naive_optimal_bbwt: rotation={naive_result}, time={naive_duration:?}");
    println!("optimal_bbwt_through_w: rotation={rotation_through_w}, time={through_w_duration:?}");
    println!("optimal_bbwt_reuse: rotation={rotation_reuse}, time={reuse_duration:?}");
    println!("optimal_bbwt_online: rotation={rotation_online}, time={online_duration:?}");

    // rotate
    let naive_rotation = [aperiodic_least_rotation[naive_result..].as_ref(), &aperiodic_least_rotation[..naive_result]].concat();
    let rotation_through_w = [aperiodic_least_rotation[rotation_through_w..].as_ref(), &aperiodic_least_rotation[..rotation_through_w]].concat();
    let rotation_reuse = [aperiodic_least_rotation[rotation_reuse..].as_ref(), &aperiodic_least_rotation[..rotation_reuse]].concat();
    let rotation_online = [aperiodic_least_rotation[rotation_online..].as_ref(), &aperiodic_least_rotation[..rotation_online]].concat();
    
    // count runs
    let runs = |text: &[u8]| bbwt(text).windows(2).filter(|w| w[0] != w[1]).count() + 1;
    let naive_runs = runs(&naive_rotation);
    assert_eq!(runs(&rotation_through_w), naive_runs);
    assert_eq!(runs(&rotation_online), naive_runs);
    assert_eq!(runs(&rotation_reuse), naive_runs);

return Ok(());

    let text_size = 7;

    // for all strings on alphabet of size text_size and length of text_size, check that the optimal bbwt implementations agree
    let mut text = vec![0; text_size];
    for i in 0..((text_size as i32).pow(text_size as u32)) {
        let mut n = i;
        for j in 0..text_size {
            text[j] = (n % text_size as i32) as u8;
            n /= text_size as i32;
        }
        if i % 10000 == 0 {
            println!("checking string {text:?} ({i})");
        }
        let (start, period_length) = shiloach(&text);
        let aperiodic_least_rotation: Vec<u8> = (0..period_length)
            .map(|k| text[(start + k) % text.len()])
            .collect();
        
        

        // compare optimal bbwt implementations
        let naive_result = optimalbbwt::naive_optimal_bbwt(&aperiodic_least_rotation);
        let rotation_through_w = optimalbbwt::optimal_bbwt_through_w(&aperiodic_least_rotation);
        let rotation_reuse = optimalbbwt::optimal_bbwt_reuse(&aperiodic_least_rotation);
        let rotation_online = optimalbbwt::online_bbwt_check(&aperiodic_least_rotation);

        // rotate
        let naive_rotation = [aperiodic_least_rotation[naive_result..].as_ref(), &aperiodic_least_rotation[..naive_result]].concat();
        let rotation_through_w = [aperiodic_least_rotation[rotation_through_w..].as_ref(), &aperiodic_least_rotation[..rotation_through_w]].concat();
        let rotation_reuse = [aperiodic_least_rotation[rotation_reuse..].as_ref(), &aperiodic_least_rotation[..rotation_reuse]].concat();
        let rotation_online = [aperiodic_least_rotation[rotation_online..].as_ref(), &aperiodic_least_rotation[..rotation_online]].concat();

        // count runs
        let runs = |text: &[u8]| bbwt(text).windows(2).filter(|w| w[0] != w[1]).count() + 1;
        let naive_runs = runs(&naive_rotation);
        assert_eq!(runs(&rotation_through_w), naive_runs);
        assert_eq!(runs(&rotation_online), naive_runs);
        assert_eq!(runs(&rotation_reuse), naive_runs);

    }

    return Ok(());

    // let mut optimal_sa = vec![0; aperiodic_least_rotation.len()];
    // optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut optimal_sa);
    // assert_eq!(
    //     optimal_sa,
    //     sais(&aperiodic_least_rotation),
    //     "suffix-array implementations disagree"
    // );

    // report_heuristic(
    //     "most repeated factors",
    //     &aperiodic_least_rotation,
    //     |input, lyn, lyn_s| most_repeated_factors(input, lyn, lyn_s),
    // );
    // report_heuristic(
    //     "most factors",
    //     &aperiodic_least_rotation,
    //     |_input, lyn, lyn_s| most_factors(lyn, lyn_s),
    // );
    // report_heuristic(
    //     "proportional repeated factors",
    //     &aperiodic_least_rotation,
    //     |input, lyn, lyn_s| most_repeated_factors_proportional(input, lyn, lyn_s),
    // );

    // report_transform("BBWT", bbwt(&text));
    // report_transform("BWT", opt_bbwt::bwt::bwt(&text));
    // Ok(())
}

fn report_heuristic(name: &str, input: &[u8], choose: impl Fn(&[u8], &[usize], &[usize]) -> usize) {
    let started = std::time::Instant::now();
    let sa = sais(input);
    let rank = sa_to_rank(&sa);
    let lyn = lyndon_table(input, &rank);
    let lyn_s = lyndon_suffix_table(input);
    let rotation = choose(input, &lyn, &lyn_s);
    let transformed = bbwt(&rotate_string(input, rotation));
    let runs = transformed
        .windows(2)
        .filter(|window| window[0] != window[1])
        .count();
    println!(
        "{name}: rotation={rotation}, runs={runs}, time={:?}",
        started.elapsed()
    );
}

fn report_transform(name: &str, transformed: Vec<u8>) {
    let runs = transformed
        .windows(2)
        .filter(|window| window[0] != window[1])
        .count();
    println!("{name}: runs={runs}");
}
