`opt-bbwt` is a modular Rust implementation of:

- the Burrows-Wheeler Transform (BWT);
- the Bijective Burrows-Wheeler Transform (BBWT);
- exact algorithms for selecting a conjugate whose BBWT has the minimum number of runs;
- heuristics for selecting a conjugate whose BBWT often has close to the minimum number of runs; and
- related algorithms including inverse transforms, suffix arrays, Lyndon factorization, Lyndon tables, Duval, and Shiloach.

This project grew out of my bachelor's thesis, [Bijective BWT: Implementing Linear-Time Lyndon Factorization of Conjugates and Discussing Open Problems](https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates). It focuses especially on the problem of selecting the conjugate of a word that yields the least number of runs in its BBWT.

Several implementations are provided for the same task to encourage experimentation and comparison, while default implementations are provided for convenience and performance in most cases.

Being a work in progress, the interface is still subject to change.

## Corpora and benchmarks

The CLI and benchmarks use Pizza&Chili corpora. Missing corpora are downloaded
on demand and cached in the workspace-level `./corpora` directory. The
`OPT_BBWT_CORPORA_DIR` environment variable can override this location.
Only the full corpus is downloaded and stored. The six supported corpus names
are `english`, `sources`, `dna`, `pitches`, `proteins`, and `dblp.xml`.

The CLI accepts either a corpus or a local input file:

```bash
cargo run -p opt-bbwt-cli -- --action bwt --corpus dna --length 50000
cargo run -p opt-bbwt-cli -- --action bwt --file ./input.bin --length 50000
```

The `sources` corpus is downloaded from Pizza&Chili's `code/sources` path.

## Documentation

The public API documentation is generated from Rust doc comments. Build and
open it locally with:

```bash
cargo doc --workspace --no-deps --open
```

To check documentation examples without opening a browser, run:

```bash
cargo test --workspace --doc
```

When publishing the library to crates.io, docs.rs automatically builds the
API documentation from these same comments.

## Optimal BBWT and conjugate selection

The project implements exact algorithms for finding a conjugate whose BBWT has the minimum number of runs. It also implements faster heuristics that often find a conjugate with a similar number of runs. In this project, an **optimal BBWT** is the BBWT obtained from a minimum-run conjugate.

The conjugate yielding the least number of runs in its BBWT is computed in a non-trivial quadratic algorithm, which has about 2x speedup over the naive brute-force approach but uses a quadratic amount of peak memory.
The heuristics choose a conjugate in O(n) time, and yield a BBWT with a number of runs that is often close to the optimal one.

The heuristics are based on the linear-time algorithm for computing the Lyndon factors of all conjugates proposed in [Bijective BWT Based Compression Schemes](https://doi.org/10.1007/s00224-025-10235-w). I first implemented this approach in C++ in [factors-of-rotations](https://github.com/TheFabbest/factors-of-rotations), then translated it to Rust and integrated it here.

## Contributing
Contributions and requests for new features or improvements are welcome! Please feel free to open an issue or submit a pull request. If you have any questions or need help, don't hesitate to reach out.

## References
- [Bijective BWT Based Compression Schemes](https://doi.org/10.1007/s00224-025-10235-w)
- [Bijective BWT: Implementing Linear-Time Lyndon Factorization of Conjugates and Discussing Open Problems](https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates)
- [factors-of-rotations](https://github.com/TheFabbest/factors-of-rotations)