This library is a modular rust implementation of the Burrows-Wheeler Transform (BWT), its bijective variant (BBWT), optimized BBWT through conjugate selection through heuristics or exact algorithms (see [my thesis](https://github.com/TheFabbest/Lyndon-Factorization-of-Conjugates)), and related algorithms (inverse BWT, inverse BBWT, duval, shiloach, lyndon tables...).

The objective of this library is to provide a fast, modular and memory-efficient implementation of the main algorithms related to the BWT. Several algorithms are implemented even for the same task, and the user can either choose the one that best fits their needs, or use the default one which is usually the fastest and most memory-efficient.

The conjugate yielding the least number of runs in its BBWT is computed in a non-trivial O(n^2) algorithm, which has about 2x speedup over the naive brute-force approach but uses a quadratic amount of peak memory.
The heuristics choose a conjugate in O(n) time, and yield a BBWT with a number of runs that is often close to the optimal one.
