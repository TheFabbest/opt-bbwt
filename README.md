This library is a modular rust implementation of:
- the Burrows-Wheeler Transform (BWT)
- the Bijective Burrows-Wheeler Transform (BBWT)
- optimized BBWT through conjugate selection by heuristics or exact algorithms (see my bachelor thesis [Bijective BWT: Implementing Linear-Time Lyndon Factorization of Conjugates and Discussing Open Problems](https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates) and the section on conjugate selection in this README)
- and related algorithms (inverse BWT, inverse BBWT, SA-IS, [optimal in-place suffix sorting](https://doi.org/10.1016/j.ic.2021.104818), Duval, Shiloach, lyndon tables...).

The objective is to provide a fast, modular and memory-efficient implementation of the main algorithms related to the BWT, this enables fast and unbiased development and evaluation of new algorithms and heuristics, while also supporting the development of new applications.
For this reason, several algorithms are implemented even for the same task, and the user can either choose the one that best fits their needs, or use the default one which is usually the fastest and most memory-efficient.

The transform APIs are grouped under `opt_bbwt::transform`: use
`transform::bwt` and `transform::bbwt` for the default implementations, and
`transform::bbwt::naive` for the reference implementation.

Being a work in progress, the interface is still subject to change, and the library is not yet to be considered stable.

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

## Conjugate selection
An interesting open problem in the context of the BBWT is the selection of a conjugate that yields the least number of runs in its BBWT. This is important because the number of runs in the BBWT is directly related to the compressibility of the transformed string, and thus to the efficiency of compression algorithms that rely on it. While our experimentation showed every rotation of the same string often yields a similar number of runs, the problem is still interesting from a theoretical point of view, and it is also relevant in practice for some specific families of strings.

The conjugate yielding the least number of runs in its BBWT is computed in a non-trivial O(n^2) algorithm, which has about 2x speedup over the naive brute-force approach but uses a quadratic amount of peak memory.
The heuristics choose a conjugate in O(n) time, and yield a BBWT with a number of runs that is often close to the optimal one.

However, we implement heuristics that often identify a good conjugate in O(n) time, and yield a BBWT with a number of runs that is often close to the optimal one. The heuristics are based on the linear-time algorithm for computing the Lyndon factors of all conjugates of a string proposed in [Bijective BWT Based Compression Schemes](https://doi.org/10.1007/s00224-025-10235-w), which was the object of my thesis and was for the first time implemented by me in C++ in [this library](https://github.com/TheFabbest/factors-of-rotations) before being translated to rust and integrated here. More precisely, our implementation of the algorithm from the paper bypasses the computation of the Lyndon trees and only relies on the Lyndon tables, which are more efficient to store, compute and use.

## Contributing
Contributions and requests for new features or improvements are welcome! Please feel free to open an issue or submit a pull request. If you have any questions or need help, don't hesitate to reach out.

## References
- [Bijective BWT Based Compression Schemes](https://doi.org/10.1007/s00224-025-10235-w)
- [Bijective BWT: Implementing Linear-Time Lyndon Factorization of Conjugates and Discussing Open Problems](https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates)
- [factors-of-rotations](https://github.com/TheFabbest/factors-of-rotations)
- [Optimal in-place suffix sorting](https://doi.org/10.1016/j.ic.2021.104818)