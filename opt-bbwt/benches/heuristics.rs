use criterion::{black_box, criterion_group, criterion_main, Criterion};
use opt_bbwt::{heuristics::{most_factors, most_repeated_factors, most_repeated_factors_proportional}, lyndon::{self, lyndon_suffix_table, lyndon_table}, rank::sa_to_rank};
use opt_bbwt::suffix_array::suffix_array;

fn heuristic_benchmarks(c: &mut Criterion) {
    let corpus_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../corpora/dna");
    let input = opt_bbwt::file::load_file_prefix(
        corpus_path.to_str().expect("corpus path must be valid UTF-8"),
        500_000,
    )
    .expect("the benchmark corpus must be available");

    let (aperiodic_least_rotation_index, period_length) = lyndon::shiloach(&input);
    let aperiodic_least_rotation = (0..period_length).map(|k| input[(aperiodic_least_rotation_index + k) % input.len()]).collect::<Vec<u8>>();
    let sa = suffix_array(&aperiodic_least_rotation);
    let rank = sa_to_rank(&sa);
    let lyn = lyndon_table(&aperiodic_least_rotation, &rank);
    let lyn_s = lyndon_suffix_table(&aperiodic_least_rotation);


    c.bench_function("most_factors", |b| b.iter(|| most_factors(black_box(&lyn), black_box(&lyn_s))));
    c.bench_function("most_repeated_factors", |b| {
        b.iter(|| most_repeated_factors(black_box(&aperiodic_least_rotation), black_box(&lyn), black_box(&lyn_s)))
    });
    c.bench_function("most_repeated_factors_proportional", |b| {
        b.iter(|| most_repeated_factors_proportional(black_box(&aperiodic_least_rotation), black_box(&lyn), black_box(&lyn_s)))
    });
}

criterion_group!(benches, heuristic_benchmarks);
criterion_main!(benches);