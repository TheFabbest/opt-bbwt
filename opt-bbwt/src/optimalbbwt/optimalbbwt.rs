use crate::bbwt::{bbwt}; use crate::lyndon::conjugatesfactors::{get_factor_with_conjugates_grouped};
use crate::lyndon::shiloach;
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
// Minimum-runs conjugate of the BBWT by a single sweep of left rotations, with certified
// in-place updates.
//
// Model (verified exhaustively in the research scripts):
//   * runs(BBWT(x)) = runs of the predecessor letters of the omega-sorted conjugates of the
//     DISTINCT Lyndon factors of x (a repeated factor only repeats letters).
//   * An "entry" is identified by a position id in the base Lyndon word s (id = abs mod n) that
//     lies in the OWNER occurrence of its factor (the occurrence whose positions are ids).
//   * left(x) = x[1..] + x[0]: the first factor v is split (Duval on v[1..]) and the moved
//     character is appended (Lyndon merging). Everything else is untouched.
//
// Certificate (the reason most entries need no work):
//   Let L = the largest LCP between two ADJACENT entries of the (physical) array. For any two
//   entries X, Y: lcp(X,Y) = min of adjacent LCPs between them <= L.
//   An entry moved by the rotation keeps its first  min(rem_old, rem_new)  characters, where rem_* is the
//   distance from its position to the end of its old / new factor occurrence (the text is the same
//   string of s). If that stable prefix is longer than L, the entry's order against every other
//   such entry (and every untouched entry) is unchanged, because their first mismatch lies inside
//   the stable prefix. So only entries within L of a factor end, or sitting on a factor start (their
//   predecessor letter changes), can move or change letter. All others are not even touched.
//   Those few entries are deleted and re-inserted by binary search; entries of vanished factors
//   that land in an occurrence which is not the owner of its new factor are deleted, and positions
//   of new owner occurrences that were not ids are inserted.
//
// Skipping with a certified lower bound (prune = true):
//   The physical array after the deletions is a SUBSEQUENCE of the true letter sequence, and runs of
//   a subsequence never exceed runs of the sequence. If runs(skeleton) >= current best, this
//   conjugate cannot beat the best, so its re-insertions are deferred (kept in a pending set). They
//   are flushed as soon as the lower bound drops below the best. Every skipped conjugate is checked
//   against the oracle in verify mode.
// =============================================================================================

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};

const MODH: u64 = (1u64 << 61) - 1;
const HBASE: u64 = 1_000_003_117_927_911;

fn mulmod(a: u64, b: u64) -> u64 {
    let t = (a as u128) * (b as u128);
    let lo = (t & MODH as u128) as u64;
    let hi = (t >> 61) as u64;
    let mut r = lo + hi;
    if r >= MODH {
        r -= MODH;
    }
    r
}

#[derive(Clone, Copy)]
struct Word {
    fs: usize,  // absolute start of the factor occurrence
    fl: usize,  // factor length
    off: usize, // offset of the conjugate inside the factor
}

#[derive(Default, Clone, Debug)]
pub struct ChainStats {
    pub steps: u64,
    pub cmp_calls: u64,
    pub cmp_chars: u64,
    pub phys_del: u64,
    pub phys_ins: u64,
    pub deferred_steps: u64,
    pub flushes: u64,
    pub rebuilds: u64,
    pub risky: u64,
    pub dead: u64,
    pub ins: u64,
    pub lmax_sum: u64,
}

#[derive(Clone, Copy)]
struct Info {
    count: usize,
    owner: usize,
}

fn cmp_words(s: &[u8], n: usize, a: &Word, b: &Word, st: &mut ChainStats) -> (usize, Ordering) {
    st.cmp_calls += 1;
    let limit = a.fl + b.fl;
    let (mut ia, mut ib) = (a.off, b.off);
    let (ba, bb) = (a.fs % n, b.fs % n);
    let mut t = 0;
    while t < limit {
        let mut pa = ba + ia;
        if pa >= n { pa -= n; }
        let mut pb = bb + ib;
        if pb >= n { pb -= n; }
        let (ca, cb) = (s[pa], s[pb]);
        if ca != cb {
            st.cmp_chars += (t + 1) as u64;
            return (t, ca.cmp(&cb));
        }
        ia += 1;
        if ia == a.fl { ia = 0; }
        ib += 1;
        if ib == b.fl { ib = 0; }
        t += 1;
    }
    st.cmp_chars += limit as u64;
    (limit, Ordering::Equal)
}

fn word_of(n: usize, id: u32, occ: &VecDeque<(usize, usize)>, ws: usize) -> Word {
    let r = (id as usize + n - ws % n) % n;
    let abs = ws + r;
    let idx = occ.partition_point(|&(st, _)| st <= abs) - 1;
    let (fs, fl) = occ[idx];
    debug_assert!(abs < fs + fl);
    Word { fs, fl, off: abs - fs }
}

fn lower_bound(
    s: &[u8], n: usize, ids: &[u32], occ: &VecDeque<(usize, usize)>, ws: usize, w: &Word, st: &mut ChainStats,
) -> usize {
    let (mut lo, mut hi) = (0usize, ids.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        let wm = word_of(n, ids[mid], occ, ws);
        if cmp_words(s, n, &wm, w, st).1 == Ordering::Less { lo = mid + 1; } else { hi = mid; }
    }
    lo
}

/// Slot of the entry `id` (whose word is `w`) in the array sorted under `occ`; never compares a word with itself.
fn locate(
    s: &[u8], n: usize, ids: &[u32], occ: &VecDeque<(usize, usize)>, ws: usize, id: u32, w: &Word, st: &mut ChainStats,
) -> usize {
    let (mut lo, mut hi) = (0usize, ids.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if ids[mid] == id { return mid; }
        let wm = word_of(n, ids[mid], occ, ws);
        if cmp_words(s, n, &wm, w, st).1 == Ordering::Less { lo = mid + 1; } else { hi = mid; }
    }
    lo
}

fn lex_less(s: &[u8], n: usize, p: (usize, usize), q: (usize, usize)) -> bool {
    let m = p.1.min(q.1);
    for k in 0..m {
        let (x, y) = (s[(p.0 + k) % n], s[(q.0 + k) % n]);
        if x != y { return x < y; }
    }
    p.1 < q.1
}

// ---- small interval-set helpers (half-open, sorted, disjoint) ----
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
    for &(al, ah) in a {
        for &(bl, bh) in b {
            let (lo, hi) = (al.max(bl), ah.min(bh));
            if lo < hi { out.push((lo, hi)); }
        }
    }
    out.sort();
    out
}
fn contains(iv: &[(usize, usize)], p: usize) -> bool {
    let k = iv.partition_point(|&(lo, _)| lo <= p);
    k > 0 && p < iv[k - 1].1
}

#[derive(Clone)]
struct Chain<'a> {
    s: &'a [u8],
    n: usize,
    i: usize, // current rotation (window start, absolute)
    pre: Vec<u64>,
    pw: Vec<u64>,
    occ: VecDeque<(usize, usize)>, // factor occurrences (abs start, len), sorted, tiling [i, i+n)
    map: HashMap<(u64, u32), Info>,
    ids: Vec<u32>,
    lt: Vec<u8>,
    lcp: Vec<u32>,
    hist: Vec<u32>,
    maxl: usize,
    diffs: usize,
    pending: Vec<u32>,
    is_pending: Vec<bool>,
    npend: usize,
    st: ChainStats,
}

impl<'a> Chain<'a> {
    fn new(s: &'a [u8]) -> Self {
        let n = s.len();
        let total = 3 * n + 2;
        let mut pre = vec![0u64; total + 1];
        let mut pw = vec![1u64; total + 1];
        for k in 0..total {
            pre[k + 1] = (mulmod(pre[k], HBASE) + s[k % n] as u64 + 1) % MODH;
            pw[k + 1] = mulmod(pw[k], HBASE);
        }
        let mut c = Chain {
            s, n, i: 0, pre, pw,
            occ: VecDeque::new(), map: HashMap::new(),
            ids: Vec::new(), lt: Vec::new(), lcp: Vec::new(),
            hist: vec![0; 2 * n + 2], maxl: 0, diffs: 0,
            pending: Vec::new(), is_pending: vec![false; n], npend: 0,
            st: ChainStats::default(),
        };
        c.occ.push_back((0, n));
        let key = c.hash_occ(0, n);
        c.map.insert(key, Info { count: 1, owner: 0 });
        c.rebuild();
        c
    }

    fn hash_occ(&self, a: usize, len: usize) -> (u64, u32) {
        let h = (self.pre[a + len] + MODH - mulmod(self.pre[a], self.pw[len])) % MODH;
        (h, len as u32)
    }

    fn runs(&self) -> usize {
        if self.ids.is_empty() { 0 } else { self.diffs + 1 }
    }

    fn letter(&self, w: &Word) -> u8 {
        let p = if w.off == 0 { w.fs + w.fl - 1 } else { w.fs + w.off - 1 };
        self.s[p % self.n]
    }

    fn hist_add(&mut self, v: u32) {
        let v = v as usize;
        if v >= self.hist.len() { self.hist.resize(v + 1, 0); }
        self.hist[v] += 1;
        if v > self.maxl { self.maxl = v; }
    }
    fn hist_rem(&mut self, v: u32) {
        let v = v as usize;
        self.hist[v] -= 1;
        while self.maxl > 0 && self.hist[self.maxl] == 0 { self.maxl -= 1; }
    }

    fn rebuild(&mut self) {
        let n = self.n;
        let ws = self.i;
        let mut items: Vec<(Word, u32)> = Vec::new();
        for (key, info) in self.map.iter() {
            for k in 0..key.1 as usize {
                let id = ((info.owner + k) % n) as u32;
                items.push((word_of(n, id, &self.occ, ws), id));
            }
        }
        let s = self.s;
        {
            let st = &mut self.st;
            items.sort_by(|a, b| cmp_words(s, n, &a.0, &b.0, st).1);
        }
        self.ids = items.iter().map(|x| x.1).collect();
        self.lt = items.iter().map(|x| self.letter(&x.0)).collect();
        self.lcp.clear();
        self.hist = vec![0; 2 * n + 2];
        self.maxl = 0;
        for j in 0..items.len().saturating_sub(1) {
            let l = cmp_words(s, n, &items[j].0, &items[j + 1].0, &mut self.st).0 as u32;
            self.lcp.push(l);
            self.hist_add(l);
        }
        self.diffs = self.lt.windows(2).filter(|w| w[0] != w[1]).count();
        for &id in &self.pending { self.is_pending[id as usize] = false; }
        self.pending.clear();
        self.npend = 0;
    }

    fn remove_slot(&mut self, j: usize) {
        let len = self.ids.len();
        let c = self.lt[j];
        if j > 0 && self.lt[j - 1] != c { self.diffs -= 1; }
        if j + 1 < len && c != self.lt[j + 1] { self.diffs -= 1; }
        if j > 0 && j + 1 < len && self.lt[j - 1] != self.lt[j + 1] { self.diffs += 1; }
        if len > 1 {
            if j == 0 {
                let v = self.lcp.remove(0);
                self.hist_rem(v);
            } else if j == len - 1 {
                let v = self.lcp.pop().unwrap();
                self.hist_rem(v);
            } else {
                let (a, b) = (self.lcp[j - 1], self.lcp[j]);
                self.hist_rem(a);
                self.hist_rem(b);
                let m = a.min(b);
                self.lcp.remove(j);
                self.lcp[j - 1] = m;
                self.hist_add(m);
            }
        }
        self.ids.remove(j);
        self.lt.remove(j);
    }

    fn insert_id(&mut self, id: u32) {
        let n = self.n;
        let ws = self.i;
        let s = self.s;
        let w = word_of(n, id, &self.occ, ws);
        let j = lower_bound(s, n, &self.ids, &self.occ, ws, &w, &mut self.st);
        let c = self.letter(&w);
        let len = self.ids.len();
        let l_left = if j > 0 {
            let wp = word_of(n, self.ids[j - 1], &self.occ, ws);
            cmp_words(s, n, &wp, &w, &mut self.st).0 as u32
        } else { 0 };
        let l_right = if j < len {
            let wn = word_of(n, self.ids[j], &self.occ, ws);
            cmp_words(s, n, &w, &wn, &mut self.st).0 as u32
        } else { 0 };
        if j > 0 && self.lt[j - 1] != c { self.diffs += 1; }
        if j < len && c != self.lt[j] { self.diffs += 1; }
        if j > 0 && j < len && self.lt[j - 1] != self.lt[j] { self.diffs -= 1; }
        if len > 0 {
            if j == 0 {
                self.lcp.insert(0, l_right);
                self.hist_add(l_right);
            } else if j == len {
                self.lcp.push(l_left);
                self.hist_add(l_left);
            } else {
                let old = self.lcp[j - 1];
                self.hist_rem(old);
                self.lcp[j - 1] = l_left;
                self.hist_add(l_left);
                self.lcp.insert(j, l_right);
                self.hist_add(l_right);
            }
        }
        self.ids.insert(j, id);
        self.lt.insert(j, c);
    }

    fn flag_pending(&mut self, id: usize) {
        if !self.is_pending[id] {
            self.is_pending[id] = true;
            self.npend += 1;
            self.pending.push(id as u32);
        }
    }

    fn flush(&mut self) {
        if self.npend == 0 { return; }
        self.st.flushes += 1;
        let list = std::mem::take(&mut self.pending);
        for id in list {
            if !self.is_pending[id as usize] { continue; }
            self.is_pending[id as usize] = false;
            self.insert_id(id);
            self.st.phys_ins += 1;
        }
        self.npend = 0;
    }

    fn duval_abs(&self, a: usize, len: usize) -> Vec<(usize, usize)> {
        let g = |k: usize| self.s[(a + k) % self.n];
        let mut out = Vec::new();
        let mut i = 0;
        while i < len {
            let mut j = i + 1;
            let mut k = i;
            while j < len && g(k) <= g(j) {
                if g(k) < g(j) { k = i; } else { k += 1; }
                j += 1;
            }
            let l = j - k;
            while i <= k {
                out.push((a + i, l));
                i += l;
            }
        }
        out
    }

    /// One left rotation: window start i -> i+1.
    fn step(&mut self, prune: bool, best: usize) {
        let n = self.n;
        let i = self.i;
        self.st.steps += 1;
        let old_occ = self.occ.clone();

        // ---- 1. factor structure of the rotated string ----
        let mut removed: Vec<(usize, usize)> = Vec::new();
        let mut added: Vec<(usize, usize)> = Vec::new();
        let (fa, fl) = self.occ.pop_front().unwrap();
        debug_assert_eq!(fa, i);
        removed.push((fa, fl));
        if fl > 1 {
            let ps = self.duval_abs(fa + 1, fl - 1);
            for &p in ps.iter().rev() { self.occ.push_front(p); }
            added.extend(ps.iter().copied());
        }
        self.occ.push_back((i + n, 1));
        added.push((i + n, 1));
        while self.occ.len() >= 2 {
            let l = self.occ[self.occ.len() - 1];
            let p = self.occ[self.occ.len() - 2];
            if lex_less(self.s, n, p, l) {
                self.occ.pop_back();
                self.occ.pop_back();
                for q in [l, p] {
                    if let Some(pos) = added.iter().position(|&x| x == q) { added.swap_remove(pos); } else { removed.push(q); }
                }
                let m = (p.0, p.1 + l.1);
                self.occ.push_back(m);
                added.push(m);
            } else { break; }
        }

        // ---- 2. distinct-factor bookkeeping ----
        let mut delta: HashMap<(u64, u32), (Vec<usize>, Vec<usize>)> = HashMap::new();
        for &(a, l) in &removed { let k = self.hash_occ(a, l); delta.entry(k).or_default().0.push(a); }
        for &(a, l) in &added { let k = self.hash_occ(a, l); delta.entry(k).or_default().1.push(a); }
        let mut vanished: Vec<(usize, usize)> = Vec::new();
        let mut appeared: Vec<(usize, usize)> = Vec::new();
        let mut need_rebuild = false;
        for (key, (rm, ad)) in delta.iter() {
            let before = self.map.get(key).map(|x| x.count).unwrap_or(0);
            let after = before + ad.len() - rm.len();
            if before == 0 {
                let owner = *ad.iter().max().unwrap();
                self.map.insert(*key, Info { count: after, owner });
                appeared.push((owner, key.1 as usize));
            } else if after == 0 {
                let info = self.map.remove(key).unwrap();
                if !rm.contains(&info.owner) { need_rebuild = true; }
                vanished.push((info.owner, key.1 as usize));
            } else {
                let info = self.map.get_mut(key).unwrap();
                if rm.contains(&info.owner) { need_rebuild = true; }
                info.count = after;
            }
        }
        self.i = i + 1;
        if need_rebuild {
            self.st.rebuilds += 1;
            self.rebuild();
            return;
        }

        // ---- 3. id-set difference in the new window [i+1, i+n+1) ----
        let nrm = |p: usize| if p == i { i + n } else { p };
        let mut vint: Vec<(usize, usize)> = Vec::new();
        for &(a, len) in &vanished {
            if a == i {
                vint.push((i + n, i + n + 1));
                if len > 1 { vint.push((i + 1, i + len)); }
            } else { vint.push((a, a + len)); }
        }
        let mut wint: Vec<(usize, usize)> = appeared.iter().map(|&(a, len)| (a, a + len)).collect();
        vint.sort();
        wint.sort();
        let dead = sub(&vint, &wint);
        let insi = sub(&wint, &vint);
        let kept = inter(&vint, &wint);

        // ---- 4. risky kept entries: stable prefix <= L, or on a factor start ----
        let lmax = self.maxl;
        self.st.lmax_sum += lmax as u64;
        let mut cand: Vec<usize> = Vec::new();
        for &(a, len) in &vanished {
            let lo = (a + len).saturating_sub(lmax).max(a);
            for p in lo..a + len { cand.push(nrm(p)); }
            cand.push(nrm(a));
        }
        for &(a, len) in &appeared {
            let lo = (a + len).saturating_sub(lmax).max(a);
            for p in lo..a + len { cand.push(p); }
            cand.push(a);
        }
        cand.sort_unstable();
        cand.dedup();
        let risky: Vec<usize> = cand.into_iter().filter(|&p| contains(&kept, p)).collect();
        self.st.risky += risky.len() as u64;

        // ---- 5. deletions under the OLD structure ----
        let mut to_delete: Vec<(usize, bool)> = Vec::new();
        for iv in &dead { for p in iv.0..iv.1 { to_delete.push((p % n, false)); } }
        self.st.dead += to_delete.len() as u64;
        for &p in &risky { to_delete.push((p % n, true)); }
        for &(id, re) in &to_delete {
            if self.is_pending[id] {
                if !re { self.is_pending[id] = false; self.npend -= 1; }
                continue;
            }
            let w = word_of(n, id as u32, &old_occ, i);
            let j = locate(self.s, n, &self.ids, &old_occ, i, id as u32, &w, &mut self.st);
            assert!(j < self.ids.len() && self.ids[j] == id as u32, "slot mismatch");
            self.remove_slot(j);
            self.st.phys_del += 1;
            if re { self.flag_pending(id); }
        }
        for iv in &insi {
            for p in iv.0..iv.1 { self.flag_pending(p % n); }
            self.st.ins += (iv.1 - iv.0) as u64;
        }

        // ---- 6. re-insertion under the NEW structure, unless certified hopeless ----
        if prune && self.npend > 0 && self.runs() >= best {
            self.st.deferred_steps += 1;
        } else {
            self.flush();
        }
    }
}

// ------------------------------- oracle & invariants (verification) -------------------------------

fn duval_slice(s: &[u8]) -> Vec<(usize, usize)> {
    let n = s.len();
    let mut out = Vec::new();
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        let mut k = i;
        while j < n && s[k] <= s[j] {
            if s[k] < s[j] { k = i; } else { k += 1; }
            j += 1;
        }
        let l = j - k;
        while i <= k { out.push((i, l)); i += l; }
    }
    out
}

fn omega_cmp_slices(a: &[u8], b: &[u8]) -> Ordering {
    let l = a.len() + b.len();
    a.iter().cycle().take(l).cmp(b.iter().cycle().take(l))
}

/// Independent oracle: letters of the sorted conjugates of the distinct Lyndon factors of x.
pub fn oracle_letters(x: &[u8]) -> Vec<u8> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut items: Vec<(Vec<u8>, u8)> = Vec::new();
    for (st, l) in duval_slice(x) {
        let f = x[st..st + l].to_vec();
        if !seen.insert(f.clone()) { continue; }
        for k in 0..l {
            let mut c = f[k..].to_vec();
            c.extend_from_slice(&f[..k]);
            items.push((c, f[(k + l - 1) % l]));
        }
    }
    items.sort_by(|a, b| omega_cmp_slices(&a.0, &b.0));
    items.into_iter().map(|e| e.1).collect()
}

pub fn runs_of_letters(v: &[u8]) -> usize {
    if v.is_empty() { 0 } else { 1 + v.windows(2).filter(|w| w[0] != w[1]).count() }
}

impl<'a> Chain<'a> {
    /// Full structural check of a (flushed) state.
    fn check_invariants(&mut self) -> Result<(), String> {
        let (n, ws, s) = (self.n, self.i, self.s);
        if self.ids.len() != self.lt.len() { return Err("ids/lt length".into()); }
        let mut st = ChainStats::default();
        let mut maxl = 0usize;
        for j in 0..self.ids.len() {
            let w = word_of(n, self.ids[j], &self.occ, ws);
            if self.letter(&w) != self.lt[j] { return Err(format!("letter mismatch at slot {}", j)); }
            if j + 1 < self.ids.len() {
                let w2 = word_of(n, self.ids[j + 1], &self.occ, ws);
                let (l, o) = cmp_words(s, n, &w, &w2, &mut st);
                if o != Ordering::Less { return Err(format!("order violated at slot {}", j)); }
                if self.lcp[j] as usize != l { return Err(format!("lcp wrong at slot {}: {} vs {}", j, self.lcp[j], l)); }
                maxl = maxl.max(l);
            }
        }
        if maxl != self.maxl { return Err(format!("maxl {} vs {}", self.maxl, maxl)); }
        let d = self.lt.windows(2).filter(|w| w[0] != w[1]).count();
        if d != self.diffs { return Err("diffs".into()); }
        Ok(())
    }
}

/// Sweep all conjugates s[i..]+s[..i] of the Lyndon word `s` (primitive). Returns
/// (index of the first conjugate with the fewest BBWT runs, that run count, statistics).
pub fn chain_min_runs(s: &[u8], prune: bool, verify: bool) -> (usize, usize, ChainStats) {
    let n = s.len();
    if n == 1 { return (0, 1, ChainStats::default()); }
    let mut c = Chain::new(s);
    let mut best = c.runs();
    let mut best_shift = 0;
    if verify {
        let o = oracle_letters(s);
        assert_eq!(c.lt, o, "initial state");
    }
    for _ in 0..n - 1 {
        let best_before = best;
        c.step(prune, best);
        if c.npend == 0 && c.runs() < best {
            best = c.runs();
            best_shift = c.i;
        }
        if verify {
            let x: Vec<u8> = (0..n).map(|k| s[(c.i + k) % n]).collect();
            let o = oracle_letters(&x);
            let true_runs = runs_of_letters(&o);
            let mut t = c.clone();
            t.flush();
            t.check_invariants().unwrap_or_else(|e| panic!("invariant after step to {}: {}", c.i, e));
            assert_eq!(t.lt, o, "letters differ at rotation {} (pending {}) s={}", c.i, c.npend, String::from_utf8_lossy(s));
            assert_eq!(t.runs(), true_runs);
            if c.npend > 0 {
                assert!(c.runs() <= true_runs, "skeleton runs exceed the true runs");
                assert!(true_runs >= best_before, "deferred a conjugate that would beat the best: {} < {}", true_runs, best_before);
            } else {
                assert_eq!(c.runs(), true_runs);
            }
        }
    }
    c.flush();
    (best_shift, best, c.st.clone())
}

// ------------------------------------- public entry points -------------------------------------

/// Index (in the conjugates of the least rotation of the primitive root) of the first conjugate
/// with the fewest BBWT runs. No verification; `prune` enables the certified skipping.
pub fn optimal_bbwt_chain(text: &[u8], prune: bool) -> usize {
    let (start, period) = shiloach(text);
    let s: Vec<u8> = (0..period).map(|k| text[(start + k) % text.len()]).collect();
    chain_min_runs(&s, prune, false).0
}

/// Same, but cross-checks every step against an independent from-scratch oracle (O(n^3 log n),
/// use for n up to a few hundred) and prints the statistics.
pub fn online_bbwt_check(input_s: &[u8]) -> usize {
    let (start, period) = shiloach(input_s);
    let s: Vec<u8> = (0..period).map(|k| input_s[(start + k) % input_s.len()]).collect();
    let verify = s.len() <= 400;
    let (shift, runs, st) = chain_min_runs(&s, true, verify);
    println!(
        "n={} verified={} best shift={} runs={} | per step: risky {:.1} dead {:.1} ins {:.1} deleted {:.1} reinserted {:.1} cmps {:.0} | deferred {} of {} steps, rebuilds {}",
        s.len(), verify, shift, runs,
        st.risky as f64 / st.steps.max(1) as f64, st.dead as f64 / st.steps.max(1) as f64,
        st.ins as f64 / st.steps.max(1) as f64, st.phys_del as f64 / st.steps.max(1) as f64,
        st.phys_ins as f64 / st.steps.max(1) as f64, st.cmp_calls as f64 / st.steps.max(1) as f64,
        st.deferred_steps, st.steps, st.rebuilds
    );
    shift
}