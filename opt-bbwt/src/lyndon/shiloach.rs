use std::cmp::Ordering;
use super::CanonicalRotationFinder;

pub struct Shiloach;

impl CanonicalRotationFinder for Shiloach {
    fn find<T: Ord>(&self, s: &[T]) -> (usize, usize) {
        if s.is_empty() {
            return (0, 0);
        }
        let (s, alpha, beta, d) = phase_one(s);
        (phase_two(s, alpha, beta), d)
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Compare s[a1+k] with s[a2+k] cyclically for k in 0..limit.
/// Returns (comparisons made, Some(order of s[a1..] vs s[a2..])) on the first
/// inequality, or (limit, None) if `limit` equalities were found.
/// Conventions: `i_max` is both "number of excluded starting points in I" and the
/// 0-based start of the next compare (the paper's `i_max + 1`, 1-based).
fn compare_circular<T: Ord>(s: &[T], a1: usize, a2: usize, limit: usize) -> (usize, Option<Ordering>) {
    let n = s.len();
    let (mut i, mut j) = (a1 % n, a2 % n);          // one modulo per call, not per element
    for k in 0..limit {
        match s[i].cmp(&s[j]) {
            Ordering::Equal => {}
            o => return (k + 1, Some(o)),
        }
        i += 1; if i == n { i = 0; }
        j += 1; if j == n { j = 0; }
    }
    (limit, None)
}

/// Effect of an inequality found after `p` comparisons: if A_i wins, j..j+p-1 are
/// excluded, otherwise i..i+p-1.  Handles wrap-around / spill-over.
fn exclude(n: usize, mut i_max: usize, mut j_max: usize, p: usize, i_wins: bool) -> (usize, usize) {
    let fl = n / 2;
    if i_wins {
        j_max += p;
        if j_max > n { i_max = i_max.max(j_max - n); j_max = n; } // wrapped into 1..
    } else {
        i_max += p;
        if i_max >= fl { j_max = j_max.max(i_max).min(n); }       // spilled into J
    }
    (i_max, j_max)
}

/// Phase 2 range (alpha, beta) after I ⊆ EX or J ⊆ EX.
fn rest(n: usize, i_max: usize, j_max: usize) -> (usize, usize) {
    if j_max >= n { (i_max, n / 2) } else { (j_max, n) }
}

/// Handles the case where ceil(n/2) equalities were found at (i_max, j_max).
fn special<T: Ord>(s: &[T], i_max: usize, j_max: usize) -> (usize, usize, usize) {
    let n = s.len();
    let fl = n / 2;
    let ce = n - fl;
    let d = gcd(n, j_max - i_max);                    // paper's unified algorithm, not n-2k
    let jp = j_max - fl;                              // j'_max
    let s_limit = n - d.max(i_max + 1).max(jp + 1);   // plain subtraction: must never underflow

    // the first ce equalities are already known; only continue if s_limit needs more
    if s_limit > ce {
        let (p, ord) = compare_circular(s, i_max + ce, j_max + ce, s_limit - ce);
        if let Some(o) = ord {
            let (i2, j2) = exclude(n, i_max, j_max, ce + p, o == Ordering::Less);
            let (a, b) = rest(n, i2, j2);
            return (a, b, n);
        }
    }
    if s_limit == n - d {                             // q > 1, odd
        if i_max >= jp { (i_max, d, d) } else { (j_max, fl + d, d) }
    } else if i_max >= jp {                           // q = 1
        (i_max.max(j_max - i_max), fl, n)
    } else {
        (j_max.max(i_max + n - jp), n, n)
    }
}

/// Returns (possibly shortened slice, alpha, beta, d).
fn phase_one<'a, T: Ord>(mut s: &'a [T]) -> (&'a [T], usize, usize, usize) {
    loop {
        let n = s.len();
        if n == 1 { return (s, 0, 1, 1); }
        let fl = n / 2;
        let ce = n - fl;
        let (mut i_max, mut j_max) = (0, fl);
        while i_max < fl && j_max < n {
            match compare_circular(s, i_max, j_max, ce) {
                (p, Some(o)) => (i_max, j_max) = exclude(n, i_max, j_max, p, o == Ordering::Less),
                (_, None) => {
                    if n % 2 == 0 && i_max == 0 && j_max == fl {
                        s = &s[..fl];                 // A = A'A': loop on the half, O(1) space
                        break;                        // restarts the outer loop
                    }
                    let (a, b, d) = special(s, i_max, j_max);
                    return (s, a, b, d);
                }
            }
        }
        if s.len() == n {                             // not halved: normal exit
            let (a, b) = rest(n, i_max, j_max);
            return (s, a, b, n);
        }
    }
}

fn phase_two<T: Ord>(s: &[T], alpha: usize, beta: usize) -> usize {
    let mut gamma = alpha;
    let mut delta = alpha + 1;
    while delta < beta {
        let limit = beta - gamma - 1;                 // paper's (beta - gamma), 1-based
        let (p, ord) = compare_circular(s, gamma, delta, limit);
        if ord == Some(Ordering::Greater) {           // A_delta wins
            gamma = delta.max(gamma + p);
            delta = gamma + 1;
        } else {
            delta += p;
        }
    }
    gamma
}
