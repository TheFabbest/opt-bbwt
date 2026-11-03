use clap::Parser;
use opt_bbwt::{
    bbwt::bbwt, heuristics::{most_factors, most_repeated_factors, most_repeated_factors_proportional}, lyndon::{lyndon_suffix_table, lyndon_table, shiloach}, rank::sa_to_rank, saca::{optimal_suffix_array_bytes, sais}, string::rotate_string,
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
    if text.len() != aperiodic_least_rotation.len() {
        eprintln!("warning: canonical rotation length differs from corpus length");
    }

    let mut optimal_sa = vec![0; aperiodic_least_rotation.len()];
    optimal_suffix_array_bytes(&aperiodic_least_rotation, &mut optimal_sa);
    assert_eq!(
        optimal_sa,
        sais(&aperiodic_least_rotation),
        "suffix-array implementations disagree"
    );

    report_heuristic(
        "most repeated factors",
        &aperiodic_least_rotation,
        |input, lyn, lyn_s| most_repeated_factors(input, lyn, lyn_s),
    );
    report_heuristic(
        "most factors",
        &aperiodic_least_rotation,
        |_input, lyn, lyn_s| most_factors(lyn, lyn_s),
    );
    report_heuristic(
        "proportional repeated factors",
        &aperiodic_least_rotation,
        |input, lyn, lyn_s| most_repeated_factors_proportional(input, lyn, lyn_s),
    );

    report_transform("BBWT", bbwt(&text));
    report_transform("BWT", opt_bbwt::bwt::bwt(&text));
    Ok(())
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
