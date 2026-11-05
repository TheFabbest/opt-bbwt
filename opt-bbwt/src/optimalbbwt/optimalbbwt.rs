use crate::bbwt::{bbwt}; use crate::lyndon::conjugatesfactors::{get_factor_with_conjugates_grouped};
use crate::lyndon::{shiloach};
use crate::saca::circularsais::optimal_csa;
use crate::saca::{circular_suffix_array};

pub fn naive_optimal_bbwt(text: &[u8]) -> usize {
    let (start, n) = shiloach(text);

    let s = (0..n).map(|k| text[(start + k) % text.len()]).collect::<Vec<u8>>();

    let mut min_runs_index = 0;
    let mut min_runs = usize::MAX;
    for i in 0..n {
        let conjugate: Vec<u8> = s[i..].iter().chain(&s[..i]).copied().collect();
        let bbwt_result = bbwt(&conjugate);
        let runs = bbwt_result.windows(2).filter(|w| w[0] != w[1]).count() + 1; // count runs in BBWT
        if runs < min_runs {
            min_runs = runs;
            min_runs_index = i;
        }
    }

    min_runs_index
}

pub fn optimal_bbwt_through_w(text: &[u8]) -> usize {
    let (start, period_length) = shiloach(text);
    let s = &(0..period_length).map(|k| text[(start + k) % text.len()]).collect::<Vec<u8>>();

    // Directly get the flattened, sorted arena
    let grouped = get_factor_with_conjugates_grouped(s);
    let sorted = grouped.sorted_factors;

    // W = concatenation
    let mut w: Vec<u8> = Vec::new();
    let mut starts = Vec::with_capacity(sorted.len());
    for &f in &sorted {  
        starts.push(w.len());
        w.extend_from_slice(f);
    }

    let mut factor_of = vec![0usize; w.len()];
    for (id, f) in sorted.iter().enumerate() {
        for p in starts[id]..starts[id] + f.len() {
            factor_of[p] = id;
        }
    }

    // Compute the circular suffix array of W
    let circular_sa = circular_suffix_array(&w);
    
    // Use 256 as a sentinel for "None" to avoid Option<u8> memory overhead and branching
    let mut run_state: Vec<(usize, u16)> = vec![(0, 256u16); s.len()];
    
    for &p in &circular_sa {
        let id = factor_of[p];
        let f = sorted[id];
        
        // REMOVED MODULO: Calculate the predecessor index safely
        let offset = p - starts[id];
        let pred = if offset == 0 {
            f[f.len() - 1]
        } else {
            f[offset - 1]
        } as u16;
        
        let idx_start = grouped.factor_starts[id];
        let idx_len = grouped.factor_lengths[id];
        let indices = &grouped.factor_indices[idx_start .. idx_start + idx_len];
        
        for &c in indices {
            let state = &mut run_state[c];
            if state.1 != pred {
                state.0 += 1;
                state.1 = pred;
            }
        }
    }

    run_state.iter().enumerate().min_by_key(|&(_, state)| state.0).map(|(i, _)| i).unwrap()
}

pub fn optimal_bbwt_reuse(text: &[u8]) -> usize {
    optimal_csa(text)
}


// =============================================================================================
// Minimum-runs conjugate of the BBWT by one sweep of left rotations (version 2).
//
// Facts used
//   * runs(BBWT(x)) = runs of the predecessor letters of the omega-sorted conjugates of the
//     DISTINCT Lyndon factors of x.
//   * Lyndon factors are non-increasing, so equal factors are adjacent ("blocks"). The input is
//     aperiodic, so the front block and the back block are never the same factor.
//   * Rotation x -> x[1..] + x[0]:
//        front : the first factor v loses its first char and splits into pieces (all strictly
//                greater than every factor that follows, so they are always NEW factors);
//        back  : the moved char is appended and merges with whole blocks at the back; a block
//                that is touched at the back dies completely.
//   * An entry is a text position (id = absolute position mod n). The entries of a factor always
//     live in the LAST occurrence of its block: the front only ever consumes the FIRST occurrence,
//     so nothing needs to move, except when the new back factor equals the previous factor (the
//     block grows at the back): then the |g| ids are shifted by |g| ("growth", rare).
//
// Certificate (why most entries are not touched)
//   L = largest LCP of two adjacent entries; any two entries share at most L characters.
//   An entry whose word changes only after p > L characters keeps its order against everything,
//   where p = min(distance to the end of its old factor, distance to the end of its new factor).
//   So only entries within L of a factor end, or on a factor start (letter changes), are deleted
//   and re-inserted by binary search.
//
// Skipping: after the deletions the array is a subsequence of the true letter sequence and
//   runs(subsequence) <= runs(sequence). If it is already >= best, the conjugate cannot win and
//   its re-insertions stay pending until the bound drops below best.
//
// Storage: entries live in a chunked sorted list (insert/delete cost O(chunk), not O(n)); every
//   entry stores the LCP with its successor, so deleting is min(lcp_prev, lcp_self) and the
//   maximum LCP comes from a histogram.
// =============================================================================================

use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};

const MODH: u64 = (1u64 << 61) - 1;
const HBASE: u64 = 1_000_003_117_927_911;
const CHUNK: usize = 512;

fn mulmod(a: u64, b: u64) -> u64 {
    let t = (a as u128) * (b as u128);
    let mut r = (t & MODH as u128) as u64 + (t >> 61) as u64;
    if r >= MODH { r -= MODH; }
    r
}

type Occ = (usize, usize); // (absolute start, length) of a factor occurrence

/// x mod n for x < 3n without a division (this sits in the innermost comparison loop).
#[inline(always)]
fn wrap(x: usize, n: usize) -> usize {
    if x < n { x } else if x < 2 * n { x - n } else { x - 2 * n }
}

#[derive(Clone, Copy)]
struct Word { fs: usize, fl: usize, off: usize } // conjugate `off` of the factor at [fs, fs+fl)

#[derive(Default, Clone, Debug)]
pub struct ChainStats {
    pub steps: u64, pub cmp_calls: u64, pub cmp_chars: u64,
    pub phys_del: u64, pub phys_ins: u64, pub deferred_steps: u64, pub flushes: u64,
    pub risky: u64, pub dead: u64, pub ins: u64, pub growth: u64, pub lmax_sum: u64,
}

/// (LCP length, order) of two conjugates, compared as infinite periodic words.
fn cmp_words(s: &[u8], n: usize, a: &Word, b: &Word, st: &mut ChainStats) -> (usize, Ordering) {
    st.cmp_calls += 1;
    let limit = a.fl + b.fl;
    let (mut ia, mut ib) = (a.off, b.off);
    let (ba, bb) = (wrap(a.fs, n), wrap(b.fs, n));
    for t in 0..limit {
        let mut pa = ba + ia; if pa >= n { pa -= n; }
        let mut pb = bb + ib; if pb >= n { pb -= n; }
        let (ca, cb) = (s[pa], s[pb]);
        if ca != cb { st.cmp_chars += (t + 1) as u64; return (t, ca.cmp(&cb)); }
        ia += 1; if ia == a.fl { ia = 0; }
        ib += 1; if ib == b.fl { ib = 0; }
    }
    st.cmp_chars += limit as u64;
    (limit, Ordering::Equal)
}

/// Word of entry `id` in the structure `occ` whose window starts at `ws`.
fn word_of(n: usize, id: u32, occ: &VecDeque<Occ>, ws: usize) -> Word {
    word_at(ws + wrap(id as usize + n - wrap(ws, n), n), occ)
}

/// Word of the entry at absolute position `abs` in the structure `occ`.
#[inline(always)]
fn word_at(abs: usize, occ: &VecDeque<Occ>) -> Word {
    let (fs, fl) = occ[occ.partition_point(|&(st, _)| st <= abs) - 1];
    Word { fs, fl, off: abs - fs }
}

/// Word of entry `id` in the structure BEFORE this rotation (window start i), reconstructed from
/// the current structure plus the list of removed occurrences (sorted by start).
fn word_old(n: usize, i: usize, id: u32, removed: &[Occ], occ: &VecDeque<Occ>) -> Word {
    let abs = i + wrap(id as usize + n - wrap(i, n), n);
    let k = removed.partition_point(|&(a, _)| a <= abs);
    if k > 0 {
        let (a, l) = removed[k - 1];
        if abs < a + l { return Word { fs: a, fl: l, off: abs - a }; }
    }
    word_at(abs, occ) // abs lies in [i+1, i+n): same place in the old and the new window
}

fn lex_less(s: &[u8], n: usize, p: Occ, q: Occ) -> bool {
    for k in 0..p.1.min(q.1) {
        let (x, y) = (s[wrap(p.0 + k, n)], s[wrap(q.0 + k, n)]);
        if x != y { return x < y; }
    }
    p.1 < q.1
}

// ---- half-open interval sets ----
fn sub(a: &[(usize, usize)], b: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for &(lo, hi) in a {
        let mut cur = lo;
        for &(bl, bh) in b {
            if bh <= cur || bl >= hi { continue; }
            if bl > cur { out.push((cur, bl)); }
            cur = cur.max(bh);
            if cur >= hi { break; }
        }
        if cur < hi { out.push((cur, hi)); }
    }
    out
}
fn inter(a: &[(usize, usize)], b: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for &(al, ah) in a { for &(bl, bh) in b {
        let (lo, hi) = (al.max(bl), ah.min(bh));
        if lo < hi { out.push((lo, hi)); }
    } }
    out.sort();
    out
}
fn contains(iv: &[(usize, usize)], p: usize) -> bool {
    let k = iv.partition_point(|&(lo, _)| lo <= p);
    k > 0 && p < iv[k - 1].1
}

// ------------------------------------ sorted entry list ------------------------------------

#[derive(Clone, Copy)]
struct Entry { id: u32, lt: u8, lcp: u32 } // lcp with the NEXT entry (0 for the last one)
type Pos = (usize, usize);                 // (chunk, offset)

#[derive(Clone)]
struct Sorted { chunks: Vec<Vec<Entry>>, len: usize, diffs: usize, hist: Vec<u32>, maxl: usize }

impl Sorted {
    fn new(n: usize) -> Self { Sorted { chunks: vec![], len: 0, diffs: 0, hist: vec![0; 2 * n + 2], maxl: 0 } }
    fn runs(&self) -> usize { if self.len == 0 { 0 } else { self.diffs + 1 } }
    fn at(&self, p: Pos) -> Entry { self.chunks[p.0][p.1] }
    fn hist_add(&mut self, v: u32) {
        let v = v as usize;
        if v >= self.hist.len() { self.hist.resize(v + 1, 0); }
        self.hist[v] += 1;
        if v > self.maxl { self.maxl = v; }
    }
    fn hist_rem(&mut self, v: u32) {
        self.hist[v as usize] -= 1;
        while self.maxl > 0 && self.hist[self.maxl] == 0 { self.maxl -= 1; }
    }
    fn from_sorted(n: usize, items: &[(u32, u8)], lcps: &[u32]) -> Self {
        let mut s = Sorted::new(n);
        let all: Vec<Entry> = items.iter().enumerate()
            .map(|(j, &(id, lt))| Entry { id, lt, lcp: if j + 1 < items.len() { lcps[j] } else { 0 } }).collect();
        for &l in lcps { s.hist_add(l); }
        s.diffs = items.windows(2).filter(|w| w[0].1 != w[1].1).count();
        s.len = all.len();
        s.chunks = all.chunks(CHUNK).map(|c| c.to_vec()).collect();
        s
    }
    /// First position whose entry is not `Less` than the target (`cmp` compares entry to target).
    fn seek(&self, mut cmp: impl FnMut(u32) -> Ordering) -> Pos {
        let (mut lo, mut hi) = (0, self.chunks.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            if cmp(self.chunks[mid][0].id) == Ordering::Less { lo = mid + 1; } else { hi = mid; }
        }
        if lo == 0 { return (0, 0); }
        let ch = &self.chunks[lo - 1];
        let (mut a, mut b) = (1, ch.len());
        while a < b {
            let mid = (a + b) / 2;
            if cmp(ch[mid].id) == Ordering::Less { a = mid + 1; } else { b = mid; }
        }
        (lo - 1, a)
    }
    /// Neighbours of the gap before position p (p.1 may equal the chunk length).
    fn gap(&self, p: Pos) -> (Option<Pos>, Option<Pos>) {
        let prev = if p.1 > 0 { Some((p.0, p.1 - 1)) } else if p.0 > 0 { Some((p.0 - 1, self.chunks[p.0 - 1].len() - 1)) } else { None };
        let next = if p.1 < self.chunks[p.0].len() { Some(p) } else if p.0 + 1 < self.chunks.len() { Some((p.0 + 1, 0)) } else { None };
        (prev, next)
    }
    fn remove(&mut self, p: Pos) {
        let e = self.at(p);
        let prev = self.gap(p).0;
        let next = if p.1 + 1 < self.chunks[p.0].len() { Some((p.0, p.1 + 1)) } else if p.0 + 1 < self.chunks.len() { Some((p.0 + 1, 0)) } else { None };
        match (prev, next) {
            (Some(a), Some(b)) => {
                let (pe, ne) = (self.at(a), self.at(b));
                if pe.lt != e.lt { self.diffs -= 1; }
                if e.lt != ne.lt { self.diffs -= 1; }
                if pe.lt != ne.lt { self.diffs += 1; }
                let m = pe.lcp.min(e.lcp);
                self.hist_rem(pe.lcp); self.hist_rem(e.lcp); self.hist_add(m);
                self.chunks[a.0][a.1].lcp = m;
            }
            (Some(a), None) => {
                let pe = self.at(a);
                if pe.lt != e.lt { self.diffs -= 1; }
                self.hist_rem(pe.lcp);
                self.chunks[a.0][a.1].lcp = 0;
            }
            (None, Some(b)) => {
                if e.lt != self.at(b).lt { self.diffs -= 1; }
                self.hist_rem(e.lcp);
            }
            (None, None) => {}
        }
        self.chunks[p.0].remove(p.1);
        if self.chunks[p.0].is_empty() { self.chunks.remove(p.0); }
        self.len -= 1;
    }
    /// Insert before position p; l_left/l_right = LCP with the previous/next entry.
    fn insert(&mut self, p: Pos, id: u32, lt: u8, l_left: u32, l_right: u32) {
        if self.chunks.is_empty() {
            self.chunks.push(vec![Entry { id, lt, lcp: 0 }]);
            self.len = 1;
            return;
        }
        let (prev, next) = self.gap(p);
        let mut e = Entry { id, lt, lcp: 0 };
        match (prev, next) {
            (Some(a), Some(b)) => {
                let (pe, ne) = (self.at(a), self.at(b));
                if pe.lt != lt { self.diffs += 1; }
                if lt != ne.lt { self.diffs += 1; }
                if pe.lt != ne.lt { self.diffs -= 1; }
                self.hist_rem(pe.lcp); self.hist_add(l_left); self.hist_add(l_right);
                self.chunks[a.0][a.1].lcp = l_left;
                e.lcp = l_right;
            }
            (Some(a), None) => {
                if self.at(a).lt != lt { self.diffs += 1; }
                self.hist_add(l_left);
                self.chunks[a.0][a.1].lcp = l_left;
            }
            (None, Some(b)) => {
                if lt != self.at(b).lt { self.diffs += 1; }
                self.hist_add(l_right);
                e.lcp = l_right;
            }
            (None, None) => unreachable!(),
        }
        self.chunks[p.0].insert(p.1, e);
        self.len += 1;
        if self.chunks[p.0].len() > 2 * CHUNK {
            let tail = self.chunks[p.0].split_off(CHUNK);
            self.chunks.insert(p.0 + 1, tail);
        }
    }
}

// ------------------------------------------ engine ------------------------------------------

#[derive(Clone)]
struct Chain<'a> {
    s: &'a [u8],
    n: usize,
    i: usize,                 // window start: the current string is s[i..] + s[..i]
    pre: Vec<u64>, pw: Vec<u64>,
    occ: VecDeque<Occ>,       // factor occurrences, sorted, tiling [i, i+n)
    seq: Sorted,
    pending: Vec<u32>, is_pending: Vec<bool>, npend: usize,
    st: ChainStats,
}

impl<'a> Chain<'a> {
    /// `sa` = start indices of the rotations of `s` in sorted order (circular suffix array); when
    /// given, it replaces the initial comparison sort. s is Lyndon, so the order is strict.
    fn new(s: &'a [u8], sa: Option<&[usize]>) -> Self {
        let n = s.len();
        let total = 3 * n + 2;
        let (mut pre, mut pw) = (vec![0u64; total + 1], vec![1u64; total + 1]);
        for k in 0..total {
            pre[k + 1] = (mulmod(pre[k], HBASE) + s[k % n] as u64 + 1) % MODH;
            pw[k + 1] = mulmod(pw[k], HBASE);
        }
        let mut c = Chain { s, n, i: 0, pre, pw, occ: VecDeque::from(vec![(0, n)]), seq: Sorted::new(n),
            pending: vec![], is_pending: vec![false; n], npend: 0, st: ChainStats::default() };
        // initial state: s is Lyndon, one factor, entries = its n conjugates
        let items: Vec<(Word, u32)> = match sa {
            Some(sa) => sa.iter().map(|&p| (word_of(n, p as u32, &c.occ, 0), p as u32)).collect(),
            None => {
                let mut v: Vec<(Word, u32)> = (0..n as u32).map(|p| (word_of(n, p, &c.occ, 0), p)).collect();
                { let st = &mut c.st; v.sort_by(|a, b| cmp_words(s, n, &a.0, &b.0, st).1); }
                v
            }
        };
        let lcps: Vec<u32> = (0..n - 1).map(|j| cmp_words(s, n, &items[j].0, &items[j + 1].0, &mut c.st).0 as u32).collect();
        let list: Vec<(u32, u8)> = items.iter().map(|(w, id)| (*id, c.letter(w))).collect();
        c.seq = Sorted::from_sorted(n, &list, &lcps);
        c
    }

    fn hash(&self, o: Occ) -> u64 { (self.pre[o.0 + o.1] + MODH - mulmod(self.pre[o.0], self.pw[o.1])) % MODH }
    fn same(&self, a: Occ, b: Occ) -> bool { a.1 == b.1 && self.hash(a) == self.hash(b) }
    fn runs(&self) -> usize { self.seq.runs() }
    fn letter(&self, w: &Word) -> u8 {
        self.s[wrap(if w.off == 0 { w.fs + w.fl - 1 } else { w.fs + w.off - 1 }, self.n)]
    }
    fn flag_pending(&mut self, id: usize) {
        if !self.is_pending[id] { self.is_pending[id] = true; self.npend += 1; self.pending.push(id as u32); }
    }
    fn unflag_pending(&mut self, id: usize) {
        if self.is_pending[id] { self.is_pending[id] = false; self.npend -= 1; }
    }

    /// Insert entry `id` (new structure) at its place by binary search.
    fn insert_id(&mut self, id: u32) {
        let (n, ws, s) = (self.n, self.i, self.s);
        let w = word_of(n, id, &self.occ, ws);
        let pos = { let (occ, st) = (&self.occ, &mut self.st);
            self.seq.seek(|e| cmp_words(s, n, &word_of(n, e, occ, ws), &w, st).1) };
        let c = self.letter(&w);
        if self.seq.chunks.is_empty() { self.seq.insert(pos, id, c, 0, 0); return; }
        let (prev, next) = self.seq.gap(pos);
        let l_left = match prev { Some(p) => {
            let wp = word_of(n, self.seq.at(p).id, &self.occ, ws);
            cmp_words(s, n, &wp, &w, &mut self.st).0 as u32 } None => 0 };
        let l_right = match next { Some(q) => {
            let wn = word_of(n, self.seq.at(q).id, &self.occ, ws);
            cmp_words(s, n, &w, &wn, &mut self.st).0 as u32 } None => 0 };
        self.seq.insert(pos, id, c, l_left, l_right);
    }

    fn flush(&mut self) {
        if self.npend == 0 { return; }
        self.st.flushes += 1;
        for id in std::mem::take(&mut self.pending) {
            if !self.is_pending[id as usize] { continue; }
            self.is_pending[id as usize] = false;
            self.insert_id(id);
            self.st.phys_ins += 1;
        }
        self.npend = 0;
    }

    fn duval_abs(&self, a: usize, len: usize) -> Vec<Occ> {
        let g = |k: usize| self.s[wrap(a + k, self.n)];
        let (mut out, mut i) = (Vec::new(), 0);
        while i < len {
            let (mut j, mut k) = (i + 1, i);
            while j < len && g(k) <= g(j) { if g(k) < g(j) { k = i; } else { k += 1; } j += 1; }
            let l = j - k;
            while i <= k { out.push((a + i, l)); i += l; }
        }
        out
    }

    /// Slot of entry `id` in the structure before the rotation.
    fn locate_old(&mut self, i: usize, id: usize, removed: &[Occ]) -> Pos {
        let (n, s) = (self.n, self.s);
        let w = word_old(n, i, id as u32, removed, &self.occ);
        let (occ, st) = (&self.occ, &mut self.st);
        let (c, o) = self.seq.seek(|e| if e == id as u32 { Ordering::Equal }
            else { cmp_words(s, n, &word_old(n, i, e, removed, occ), &w, st).1 });
        let p = if o == self.seq.chunks[c].len() { (c + 1, 0) } else { (c, o) };
        assert!(self.seq.at(p).id == id as u32, "entry not found");
        p
    }

    /// One left rotation: window start i -> i+1.
    fn step(&mut self, prune: bool, best: usize) {
        let (n, i) = (self.n, self.i);
        self.st.steps += 1;

        // A. new factorization: split the first factor, append the moved char and merge at the back
        let front_dup = self.occ.len() > 1 && self.same(self.occ[0], self.occ[1]);
        let (fa, fl) = self.occ.pop_front().unwrap();
        let mut removed: Vec<Occ> = vec![(fa, fl)];
        let mut added: Vec<Occ> = Vec::new();
        if fl > 1 {
            let ps = self.duval_abs(fa + 1, fl - 1);
            for &p in ps.iter().rev() { self.occ.push_front(p); }
            added.extend(ps);
        }
        self.occ.push_back((i + n, 1));
        added.push((i + n, 1));
        while self.occ.len() >= 2 {
            let (l, p) = (self.occ[self.occ.len() - 1], self.occ[self.occ.len() - 2]);
            if !lex_less(self.s, n, p, l) { break; }
            self.occ.pop_back(); self.occ.pop_back();
            for q in [l, p] {
                if let Some(k) = added.iter().position(|&x| x == q) { added.swap_remove(k); } else { removed.push(q); }
            }
            let m = (p.0, p.1 + l.1);
            self.occ.push_back(m); added.push(m);
        }
        removed.sort(); added.sort();
        let front_persists = front_dup && removed.binary_search(&(fa + fl, fl)).is_err();

        // B. which distinct factors die / are born (owner = last occurrence of the block)
        let mut vanished: Vec<Occ> = Vec::new();
        let mut k = removed.len();
        while k > 0 {
            let rep = removed[k - 1];
            k -= 1;
            while k > 0 && self.same(removed[k - 1], rep) { k -= 1; }
            if !(front_persists && self.same(rep, (fa, fl))) { vanished.push(rep); }
        }
        // the back factor equals the previous (surviving) factor: the block grows
        let back = *self.occ.back().unwrap();
        let mut growth: Option<(Occ, Occ)> = None;
        if self.occ.len() >= 2 {
            let prev = self.occ[self.occ.len() - 2];
            if added.binary_search(&prev).is_err() && self.same(prev, back) {
                growth = Some((prev, back));
                added.retain(|&q| q != back);
            }
        }
        let mut appeared: Vec<Occ> = Vec::new();
        let mut k = 0;
        while k < added.len() {
            let mut j = k;
            while j + 1 < added.len() && added[j + 1].0 == added[j].0 + added[j].1 && self.same(added[j + 1], added[k]) { j += 1; }
            appeared.push(added[j]);
            k = j + 1;
        }
        self.i = i + 1;

        // C. entry positions that die / are born / stay, in the new window [i+1, i+n+1)
        let mut vint: Vec<(usize, usize)> = Vec::new();
        for &(a, len) in &vanished {
            if a == i { vint.push((i + n, i + n + 1)); if len > 1 { vint.push((i + 1, i + len)); } }
            else { vint.push((a, a + len)); }
        }
        let mut wint: Vec<(usize, usize)> = appeared.iter().map(|&(a, len)| (a, a + len)).collect();
        vint.sort(); wint.sort();
        let (dead, insi, kept) = (sub(&vint, &wint), sub(&wint, &vint), inter(&vint, &wint));

        // D. kept entries that can still move: within L of a factor end, or on a factor start
        let lmax = self.seq.maxl;
        self.st.lmax_sum += lmax as u64;
        let nrm = |p: usize| if p == i { i + n } else { p };
        let mut cand: Vec<usize> = Vec::new();
        for &(a, len) in &vanished {
            for p in (a + len).saturating_sub(lmax).max(a)..a + len { cand.push(nrm(p)); }
            cand.push(nrm(a));
        }
        for &(a, len) in &appeared {
            for p in (a + len).saturating_sub(lmax).max(a)..a + len { cand.push(p); }
            cand.push(a);
        }
        cand.sort_unstable(); cand.dedup();
        let risky: Vec<usize> = cand.into_iter().filter(|&p| contains(&kept, p)).collect();
        self.st.risky += risky.len() as u64;

        // E. delete dead + risky entries (old structure), queue re-insertions
        let mut to_delete: Vec<(usize, bool)> = Vec::new();
        for iv in &dead { for p in iv.0..iv.1 { to_delete.push((p % n, false)); } }
        self.st.dead += to_delete.len() as u64;
        for &p in &risky { to_delete.push((p % n, true)); }
        for &(id, requeue) in &to_delete {
            if self.is_pending[id] { if !requeue { self.unflag_pending(id); } continue; }
            let p = self.locate_old(i, id, &removed);
            self.seq.remove(p);
            self.st.phys_del += 1;
            if requeue { self.flag_pending(id); }
        }
        for iv in &insi { for p in iv.0..iv.1 { self.flag_pending(p % n); } self.st.ins += (iv.1 - iv.0) as u64; }
        if let Some((pv, bk)) = growth {   // block grows at the back: shift its ids by |g|
            self.st.growth += 1;
            let slots: Vec<(usize, usize, Option<Pos>)> = (0..bk.1).map(|k| {
                let (a, b) = ((pv.0 + k) % n, (bk.0 + k) % n);
                let p = if self.is_pending[a] { None } else { Some(self.locate_old(i, a, &removed)) };
                (a, b, p)
            }).collect();
            for (a, b, p) in slots {
                match p {
                    Some(p) => self.seq.chunks[p.0][p.1].id = b as u32,
                    None => { self.unflag_pending(a); self.flag_pending(b); }
                }
            }
        }

        // F. re-insert, unless the conjugate provably cannot beat `best`
        if prune && self.npend > 0 && self.runs() >= best { self.st.deferred_steps += 1; } else { self.flush(); }
    }
}

// ------------------------------- oracle & invariants (verification) -------------------------------

fn duval_slice(s: &[u8]) -> Vec<(usize, usize)> {
    let n = s.len();
    let (mut out, mut i) = (Vec::new(), 0);
    while i < n {
        let (mut j, mut k) = (i + 1, i);
        while j < n && s[k] <= s[j] { if s[k] < s[j] { k = i; } else { k += 1; } j += 1; }
        let l = j - k;
        while i <= k { out.push((i, l)); i += l; }
    }
    out
}

/// Independent oracle: letters of the sorted conjugates of the distinct Lyndon factors of x.
pub fn oracle_letters(x: &[u8]) -> Vec<u8> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut items: Vec<(Vec<u8>, u8)> = Vec::new();
    for (st, l) in duval_slice(x) {
        let f = x[st..st + l].to_vec();
        if !seen.insert(f.clone()) { continue; }
        for k in 0..l {
            let mut c = f[k..].to_vec(); c.extend_from_slice(&f[..k]);
            items.push((c, f[(k + l - 1) % l]));
        }
    }
    items.sort_by(|a, b| { let l = a.0.len() + b.0.len(); a.0.iter().cycle().take(l).cmp(b.0.iter().cycle().take(l)) });
    items.into_iter().map(|e| e.1).collect()
}

pub fn runs_of_letters(v: &[u8]) -> usize {
    if v.is_empty() { 0 } else { 1 + v.windows(2).filter(|w| w[0] != w[1]).count() }
}

impl<'a> Chain<'a> {
    fn letters(&self) -> Vec<u8> { self.seq.chunks.iter().flatten().map(|e| e.lt).collect() }
    fn check_invariants(&mut self) -> Result<(), String> {
        let (n, ws, s) = (self.n, self.i, self.s);
        let flat: Vec<Entry> = self.seq.chunks.iter().flatten().copied().collect();
        if flat.len() != self.seq.len { return Err("len".into()); }
        if self.seq.chunks.iter().any(|c| c.is_empty() || c.len() > 2 * CHUNK) { return Err("chunk size".into()); }
        let mut st = ChainStats::default();
        let mut maxl = 0;
        for j in 0..flat.len() {
            let w = word_of(n, flat[j].id, &self.occ, ws);
            if self.letter(&w) != flat[j].lt { return Err(format!("letter at {}", j)); }
            if j + 1 < flat.len() {
                let w2 = word_of(n, flat[j + 1].id, &self.occ, ws);
                let (l, o) = cmp_words(s, n, &w, &w2, &mut st);
                if o != Ordering::Less { return Err(format!("order at {}", j)); }
                if flat[j].lcp as usize != l { return Err(format!("lcp at {}", j)); }
                maxl = maxl.max(l);
            } else if flat[j].lcp != 0 { return Err("last lcp".into()); }
        }
        if maxl != self.seq.maxl { return Err(format!("maxl {} vs {}", self.seq.maxl, maxl)); }
        if flat.windows(2).filter(|w| w[0].lt != w[1].lt).count() != self.seq.diffs { return Err("diffs".into()); }
        Ok(())
    }
}

/// Sweep all conjugates s[i..]+s[..i] of the Lyndon word `s`. Returns (index of the first conjugate
/// with the fewest BBWT runs, that run count, statistics). `verify` cross-checks every step.
pub fn chain_min_runs(s: &[u8], prune: bool, verify: bool) -> (usize, usize, ChainStats) {
    chain_min_runs_sa(s, None, prune, verify)
}

/// Same, with the circular suffix array of `s` (sorted rotation starts) supplied by the caller.
pub fn chain_min_runs_sa(s: &[u8], sa: Option<&[usize]>, prune: bool, verify: bool) -> (usize, usize, ChainStats) {
    let n = s.len();
    if n == 1 { return (0, 1, ChainStats::default()); }
    let mut c = Chain::new(s, sa);
    let (mut best, mut best_shift) = (c.runs(), 0);
    if verify { assert_eq!(c.letters(), oracle_letters(s), "initial state"); }
    for _ in 0..n - 1 {
        let best_before = best;
        c.step(prune, best);
        if c.npend == 0 && c.runs() < best { best = c.runs(); best_shift = c.i; }
        if verify {
            let x: Vec<u8> = (0..n).map(|k| s[(c.i + k) % n]).collect();
            let o = oracle_letters(&x);
            let true_runs = runs_of_letters(&o);
            let mut t = c.clone();
            t.flush();
            t.check_invariants().unwrap_or_else(|e| panic!("invariant after step to {}: {}", c.i, e));
            assert_eq!(t.letters(), o, "letters differ at rotation {} (pending {})", c.i, c.npend);
            if c.npend > 0 {
                assert!(c.runs() <= true_runs, "skeleton runs exceed true runs");
                assert!(true_runs >= best_before, "deferred a conjugate that beats the best: {} < {}", true_runs, best_before);
            } else { assert_eq!(c.runs(), true_runs); }
        }
    }
    c.flush();
    (best_shift, best, c.st.clone())
}

// ------------------------------------- public entry points -------------------------------------

/// Index (among the conjugates of the least rotation of the primitive root) of the first conjugate
/// with the fewest BBWT runs. `prune` enables the certified skipping.
pub fn optimal_bbwt_chain(text: &[u8], prune: bool) -> usize {
    let (start, period) = shiloach(text);
    let s: Vec<u8> = (0..period).map(|k| text[(start + k) % text.len()]).collect();
    let sa = circular_suffix_array(&s); // initial order of the n rotations: the only suffix-array computation of the run
    chain_min_runs_sa(&s, Some(&sa), prune, false).0
}

/// Same sweep, cross-checked against the oracle at every step when n <= 400, with statistics.
pub fn online_bbwt_check(input_s: &[u8]) -> usize {
    let (start, period) = shiloach(input_s);
    let s: Vec<u8> = (0..period).map(|k| input_s[(start + k) % input_s.len()]).collect();
    let verify = s.len() <= 400;
    let sa = circular_suffix_array(&s);
    let (shift, runs, st) = chain_min_runs_sa(&s, Some(&sa), true, verify);
    let sn = st.steps.max(1) as f64;
    println!("n={} verified={} best shift={} runs={} | per step: risky {:.1} deleted {:.1} reinserted {:.1} cmps {:.0} | deferred {} of {} steps, block growths {}",
        s.len(), verify, shift, runs, st.risky as f64 / sn, st.phys_del as f64 / sn, st.phys_ins as f64 / sn,
        st.cmp_calls as f64 / sn, st.deferred_steps, st.steps, st.growth);
    shift
}