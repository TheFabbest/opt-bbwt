use criterion::{black_box, criterion_group, criterion_main, Criterion};
use opt_bbwt::optimal_bbwt::{
    chain::optimal_bbwt_chain, globalsort::optimal_bbwt_global_sort, naive::naive_optimal_bbwt,
    reuse::optimal_bbwt_reuse,
};

fn optimal_bbwt_benchmarks(c: &mut Criterion) {
    let input = opt_bbwt::file::load_corpus("dna", 1_500)
        .expect("the benchmark corpus must be available or downloadable");
    c.bench_function("naive_optimal_bbwt", |b| {
        b.iter(|| naive_optimal_bbwt(black_box(&input)))
    });
    c.bench_function("optimal_bbwt_reuse", |b| {
        b.iter(|| optimal_bbwt_reuse(black_box(&input)))
    });
    c.bench_function("optimal_bbwt_global_sort", |b| {
        b.iter(|| optimal_bbwt_global_sort(black_box(&input)))
    });
    c.bench_function("optimal_bbwt_chain", |b| {
        b.iter(|| optimal_bbwt_chain(black_box(&input), black_box(true)))
    });
}

criterion_group!(benches, optimal_bbwt_benchmarks);
criterion_main!(benches);
