use criterion::{black_box, criterion_group, criterion_main, Criterion};
use opt_bbwt::transform::{bbwt::bbwt, bwt::bwt};

fn transform_benchmarks(c: &mut Criterion) {
    let input = opt_bbwt::file::load_corpus("dna", 50_000)
        .expect("the benchmark corpus must be available or downloadable");
    c.bench_function("bwt", |b| b.iter(|| bwt(black_box(&input))));
    c.bench_function("bbwt", |b| b.iter(|| bbwt(black_box(&input))));
}

criterion_group!(benches, transform_benchmarks);
criterion_main!(benches);
