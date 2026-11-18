use crate::lyndon::shiloach;
use crate::suffix_array::circular_suffix_array;
use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};

// =============================================================================================
// Finds the conjugate of a string whose BBWT has the fewest runs, by sweeping ALL rotations
//     x = a·y  ->  y·a        (one "left rotation" per step)
// and updating a sorted list instead of recomputing it.
//
// VOCABULARY
//   text      the base string (a Lyndon word: the least rotation of its primitive root); text_len = its length.
//   window    the current conjugate is text[window_start..] followed by text[..window_start].
//             Positions are "absolute": window_start..window_start+text_len, so always < 2*text_len.
//             They index `doubled` = text·text directly (no wrapping); `wrap` maps one to a text index.
//   factor    one occurrence of a Lyndon factor inside the window, stored as (start, len).
//   block     several equal factors in a row (Lyndon factors are non-increasing, so equal ones are adjacent).
//   entry     one conjugate of one DISTINCT factor, e.g. factor "aab" has entries "aab","aba","baa".
//             It is identified by `id` = a text position inside the LAST occurrence of its block.
//   owner     the last occurrence of a block: the factor whose positions are the ids of the entries.
//   letter    the character just before the entry inside its factor (what the BBWT outputs for it).
//   runs      number of maximal blocks of equal letters in the entries sorted in omega-order.
//
// WHAT ONE ROTATION DOES   (c = text[window_start] is the moved character)
//   front  the first factor v loses its first char. v[1..] factorizes into "pieces", which are read
//          from the Lyndon table of `text`. Pieces are always new factors.
//   back   c is appended as a one-char factor and merges with the factors before it while
//          "second-to-last < last". The result is ONE new factor `back`.
//
// INVARIANTS (proved in the comments where they are used, checked by the oracle in `verify` mode)
//   I1  Every factor merged into `back` starts with a letter smaller than c (a Lyndon word of
//       length > 1 that ends with c starts with a letter < c, and factors are non-increasing).
//       Every piece and every copy of v starts with a letter >= c. So merging never reaches the
//       pieces or the copies of the first block, and a piece never equals `back`.
//   I2  The pieces are the Lyndon factors of v[1..]; they follow from the Lyndon table of `text`:
//       v is the longest Lyndon prefix of text[window_start..] (it never wraps), and the walk
//       p -> p + table[p] from window_start+1 lands exactly on the end of v.
//
// HOW ENTRIES CHANGE (owners die, owners are born)
//   dead owner  = the first factor v if it has no second copy; every merged block (its last copy).
//   born owner  = the last piece of each block of equal pieces; `back`, unless it equals the old
//                 factor before it: then that block only grows (its last copy is now `back`) and
//                 its entries are renamed (ids + |back|) instead of being deleted and inserted.
//                 Renaming keeps them in the list, which keeps the lower bound used for skipping high.
//   An id that is a dead owner position and not a born owner position is deleted for good;
//   one that is a born owner position and not a dead one is inserted; one that is both ("kept")
//   keeps its id but its word and/or letter changed, so it may need a new place in the order.
//   The first position of the window leaves at the front and re-enters at the back
//   (position p -> p + text_len); the other renaming is the growth of a block (see above).
//
// WHY MOST ENTRIES ARE NOT TOUCHED
//   max_lcp = largest LCP between two neighbours in the sorted list = the longest common prefix
//   any two entries can have. An entry whose word is unchanged for its first p > max_lcp characters
//   keeps its place in the order. So only kept entries within max_lcp characters of a factor end
//   (their word changes there) or on a factor start (their letter changes) are deleted and
//   re-inserted ("risky").
//
// WHY WORK CAN BE SKIPPED
//   After the deletions the list is a subsequence of the true letter sequence, and a subsequence
//   never has more runs. If that is already >= best, this conjugate cannot win: the re-insertions
//   stay "pending" until a later step has a smaller bound.
//
// STORAGE
//   Sorted entries live in chunks of ~CHUNK (insert/delete cost O(chunk)). Each entry stores its LCP
//   with the next entry, so deleting an entry gives min(lcp_prev, lcp_self) for the joined pair,
//   and a histogram of all LCPs gives max_lcp in O(1).
// =============================================================================================

const HASH_MOD: u64 = (1u64 << 61) - 1;
const HASH_BASE: u64 = 1_000_003_117_927_911;
const CHUNK: usize = 128;

/// (a * b) mod (2^61 - 1): 128-bit product, then fold the high bits onto the low ones.
fn mul_mod(a: u64, b: u64) -> u64 {
    let product = (a as u128) * (b as u128);
    let r = (product & HASH_MOD as u128) as u64 + (product >> 61) as u64;
    if r >= HASH_MOD {
        r - HASH_MOD
    } else {
        r
    }
}

/// One occurrence of a Lyndon factor in the window: doubled[start..start+len] (absolute positions).
#[derive(Clone, Copy)]
struct Factor {
    start: usize,
    len: usize,
}
impl Factor {
    fn end(&self) -> usize {
        self.start + self.len
    }
}

/// x mod text_len for 0 <= x < 2*text_len: an absolute window position -> the text index (an entry id).
#[inline(always)]
fn wrap(x: usize, text_len: usize) -> usize {
    if x >= text_len {
        x - text_len
    } else {
        x
    }
}

/// A conjugate of a factor: the factor's text rotated left by `off`. It is infinite-periodic
/// when compared (see `compare_conjugates`).
#[derive(Clone, Copy)]
struct Conjugate {
    start: usize,
    len: usize,
    off: usize,
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
    pub risky: u64,
    pub dead: u64,
    pub ins: u64,
    pub growth: u64,
    pub lmax_sum: u64,
}

/// Length of the longest common prefix of two slices of equal length, 8 bytes at a time.
#[inline(always)]
fn lcp(a: &[u8], b: &[u8]) -> usize {
    let mut i = 0;
    for (x, y) in a.chunks_exact(8).zip(b.chunks_exact(8)) {
        let diff = u64::from_le_bytes(x.try_into().unwrap()) ^ u64::from_le_bytes(y.try_into().unwrap());
        if diff != 0 {
            return i + diff.trailing_zeros() as usize / 8;
        }
        i += 8;
    }
    while i < a.len() && a[i] == b[i] {
        i += 1;
    }
    i
}

/// Compares two conjugates as infinite periodic words (each one repeats its own factor forever).
/// Returns (length of their longest common prefix, Less/Equal/Greater).
/// How it works: walk both words in parallel, each inside its own factor, comparing whole stretches
/// (up to the end of the nearer factor) with `lcp`, and wrapping to the factor's start when it ends.
/// Two periodic words with periods la, lb that agree on la+lb characters are equal (Fine–Wilf),
/// so `la + lb` characters are enough to decide.
fn compare_conjugates(t: &[u8], a: &Conjugate, b: &Conjugate, stats: &mut ChainStats) -> (usize, Ordering) {
    stats.cmp_calls += 1;
    let max_len = a.len + b.len;
    let (mut pa, mut pb) = (a.start + a.off, b.start + b.off);
    let (end_a, end_b) = (a.start + a.len, b.start + b.len);
    let mut k = 0;
    while k < max_len {
        let seg = (end_a - pa).min(end_b - pb).min(max_len - k);
        let l = lcp(&t[pa..pa + seg], &t[pb..pb + seg]);
        k += l;
        if l < seg {
            stats.cmp_chars += k as u64 + 1;
            return (k, t[pa + l].cmp(&t[pb + l]));
        }
        pa += seg;
        pb += seg;
        if pa == end_a {
            pa = a.start;
        }
        if pb == end_b {
            pb = b.start;
        }
    }
    stats.cmp_chars += max_len as u64;
    (max_len, Ordering::Equal)
}

/// The conjugate that entry `id` stands for in the window starting at `window_start`, whose
/// factors are `factors`. (id -> absolute position -> the factor containing it, by binary search on
/// the sorted factor starts.)
fn conjugate_of(text_len: usize, window_start: usize, id: u32, factors: &VecDeque<Factor>) -> Conjugate {
    let id = id as usize;
    let abs = if id >= window_start { id } else { id + text_len };
    let f = factors[factors.partition_point(|f| f.start <= abs) - 1];
    Conjugate { start: f.start, len: f.len, off: abs - f.start }
}

/// Is the factor f smaller than g as a plain string (a proper prefix is smaller)? Used for the
/// Lyndon merge: if f < g then f·g is a Lyndon word, so neighbouring factors f < g must be merged.
fn factor_is_smaller(t: &[u8], f: Factor, g: Factor) -> bool {
    let m = f.len.min(g.len);
    let l = lcp(&t[f.start..f.start + m], &t[g.start..g.start + m]);
    if l < m {
        t[f.start + l] < t[g.start + l]
    } else {
        f.len < g.len
    }
}

/// One entry of the sorted list.
#[derive(Clone, Copy)]
struct Entry {
    id: u32,       // text position that identifies the entry (see vocabulary)
    letter: u8,    // the BBWT letter of this entry
    lcp_next: u32, // longest common prefix with the NEXT entry in the list (0 for the last one)
}
/// Location in the list: (chunk index, index inside the chunk).
type Slot = (usize, usize);

/// All entries in omega-sorted order, kept in chunks, with the data needed to know the run count
/// and the largest neighbour LCP at any moment.
#[derive(Clone)]
struct SortedEntries {
    chunks: Vec<Vec<Entry>>,
    chunk_uid: Vec<u32>,     // chunk_uid[i] = the permanent name of chunks[i]
    chunk_index: Vec<usize>, // chunk_index[uid] = where that chunk is now (usize::MAX once it is gone)
    chunk_of: Vec<u32>, // chunk_of[id] = name of the chunk that holds entry `id` (or held it last)
    len: usize,
    letter_changes: usize, // number of adjacent pairs with different letters; runs = this + 1
    lcp_histogram: Vec<u32>, // lcp_histogram[v] = how many adjacent pairs have LCP exactly v
    max_lcp: usize,        // largest v with lcp_histogram[v] > 0
}

impl SortedEntries {
    fn new(text_len: usize) -> Self {
        SortedEntries {
            chunks: vec![],
            chunk_uid: vec![],
            chunk_index: vec![],
            chunk_of: vec![0; text_len],
            len: 0,
            letter_changes: 0,
            lcp_histogram: vec![0; 2 * text_len + 2],
            max_lcp: 0,
        }
    }
    fn runs(&self) -> usize {
        if self.len == 0 {
            0
        } else {
            self.letter_changes + 1
        }
    }
    fn at(&self, slot: Slot) -> Entry {
        self.chunks[slot.0][slot.1]
    }
    fn histogram_add(&mut self, v: u32) {
        let v = v as usize;
        if v >= self.lcp_histogram.len() {
            self.lcp_histogram.resize(v + 1, 0);
        }
        self.lcp_histogram[v] += 1;
        if v > self.max_lcp {
            self.max_lcp = v;
        }
    }
    fn histogram_remove(&mut self, v: u32) {
        self.lcp_histogram[v as usize] -= 1;
        while self.max_lcp > 0 && self.lcp_histogram[self.max_lcp] == 0 {
            self.max_lcp -= 1;
        }
    }
    /// Builds the list from entries that are already sorted. `lcps[j]` = LCP of entries j and j+1.
    fn from_sorted(text_len: usize, entries: &[(u32, u8)], lcps: &[u32]) -> Self {
        let mut list = SortedEntries::new(text_len);
        let all: Vec<Entry> = entries
            .iter()
            .enumerate()
            .map(|(j, &(id, letter))| Entry {
                id,
                letter,
                lcp_next: if j + 1 < entries.len() { lcps[j] } else { 0 },
            })
            .collect();
        for &l in lcps {
            list.histogram_add(l);
        }
        list.letter_changes = entries.windows(2).filter(|w| w[0].1 != w[1].1).count();
        list.len = all.len();
        for chunk in all.chunks(CHUNK) {
            list.insert_chunk(list.chunks.len(), chunk.to_vec());
        }
        list
    }
    /// Adds a chunk at position `index` and records, for every entry in it, which chunk it is in.
    fn insert_chunk(&mut self, index: usize, chunk: Vec<Entry>) {
        let uid = self.chunk_index.len() as u32;
        for e in &chunk {
            self.chunk_of[e.id as usize] = uid;
        }
        self.chunk_index.push(index);
        self.chunk_uid.insert(index, uid);
        self.chunks.insert(index, chunk);
        for i in index + 1..self.chunks.len() {
            self.chunk_index[self.chunk_uid[i] as usize] = i; // chunks after it moved one place
        }
    }
    fn remove_chunk(&mut self, index: usize) {
        self.chunk_index[self.chunk_uid[index] as usize] = usize::MAX;
        self.chunk_uid.remove(index);
        self.chunks.remove(index);
        for i in index..self.chunks.len() {
            self.chunk_index[self.chunk_uid[i] as usize] = i;
        }
    }
    /// The slot of entry `id`, which must be in the list. No comparisons: the chunk is looked up in
    /// `chunk_of`, then the id is searched for inside the chunk (at most 2*CHUNK ids).
    fn find(&self, id: u32) -> Slot {
        let chunk = self.chunk_index[self.chunk_of[id as usize] as usize];
        let offset = self.chunks[chunk]
            .iter()
            .position(|e| e.id == id)
            .expect("entry not in the list");
        (chunk, offset)
    }
    /// Gives the entry at `slot` a new id (same word, same place in the order).
    fn rename(&mut self, slot: Slot, new_id: u32) {
        self.chunks[slot.0][slot.1].id = new_id;
        self.chunk_of[new_id as usize] = self.chunk_uid[slot.0];
    }
    /// Binary search: first slot whose entry is NOT smaller than the target.
    /// `compare_to_target(id)` must say how the entry `id` compares to the target (Less/Equal/Greater).
    /// Two levels: first over the first entry of every chunk, then inside one chunk.
    fn seek(&self, mut compare_to_target: impl FnMut(u32) -> Ordering) -> Slot {
        let (mut lo, mut hi) = (0, self.chunks.len());
        while lo < hi {
            let mid = (lo + hi) / 2;
            if compare_to_target(self.chunks[mid][0].id) == Ordering::Less {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if lo == 0 {
            return (0, 0);
        }
        let chunk = &self.chunks[lo - 1];
        let (mut a, mut b) = (1, chunk.len());
        while a < b {
            let mid = (a + b) / 2;
            if compare_to_target(chunk[mid].id) == Ordering::Less {
                a = mid + 1;
            } else {
                b = mid;
            }
        }
        (lo - 1, a)
    }
    /// Same result as `seek`, for a target that is probably close to where the entry `id` used to be:
    /// `offset` is its old index inside the chunk it was in. The search gallops outward from there
    /// (1, 2, 4, ... steps, then a binary search inside the bracket), so it costs about
    /// log2(distance) comparisons instead of log2(list length). If the answer is not inside that
    /// chunk (or the chunk is gone) it falls back to `seek`.
    fn seek_near(
        &self,
        id: u32,
        offset: usize,
        mut compare_to_target: impl FnMut(u32) -> Ordering,
    ) -> Slot {
        let c = self.chunk_index[self.chunk_of[id as usize] as usize];
        if c == usize::MAX {
            return self.seek(compare_to_target);
        }
        let chunk = &self.chunks[c];
        let guess = offset.min(chunk.len());
        // the answer is the first index in [lo, hi] whose entry is not smaller than the target
        let (mut lo, mut hi) = (0, chunk.len());
        if guess < hi && compare_to_target(chunk[guess].id) == Ordering::Less {
            lo = guess + 1; // gallop right
            let mut step = 1;
            while lo < hi {
                let probe = (lo + step - 1).min(hi - 1);
                if compare_to_target(chunk[probe].id) == Ordering::Less {
                    lo = probe + 1;
                    step *= 2;
                } else {
                    hi = probe;
                    break;
                }
            }
        } else {
            hi = guess; // gallop left
            let mut step = 1;
            while hi > 0 {
                let probe = hi.saturating_sub(step);
                if compare_to_target(chunk[probe].id) == Ordering::Less {
                    lo = probe + 1;
                    break;
                }
                hi = probe;
                step *= 2;
            }
        }
        while lo < hi {
            let mid = (lo + hi) / 2;
            if compare_to_target(chunk[mid].id) == Ordering::Less {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        // the bracket may have run into a chunk border: then the neighbour chunk decides
        let left_ok = lo > 0
            || c == 0
            || compare_to_target(self.chunks[c - 1].last().unwrap().id) == Ordering::Less;
        let right_ok = lo < chunk.len()
            || c + 1 == self.chunks.len()
            || compare_to_target(self.chunks[c + 1][0].id) != Ordering::Less;
        if left_ok && right_ok {
            (c, lo)
        } else {
            self.seek(compare_to_target)
        }
    }
    /// The first slot at or after `slot` (None past the end); `slot.1` may equal the chunk length.
    fn first_from(&self, slot: Slot) -> Option<Slot> {
        if slot.1 < self.chunks.get(slot.0).map_or(0, Vec::len) {
            Some(slot)
        } else if slot.0 + 1 < self.chunks.len() {
            Some((slot.0 + 1, 0))
        } else {
            None
        }
    }
    /// The slots just before and just after the gap in front of `slot` (None at the list ends).
    fn neighbours_of_gap(&self, slot: Slot) -> (Option<Slot>, Option<Slot>) {
        let prev = match slot {
            (c, i) if i > 0 => Some((c, i - 1)),
            (c, _) if c > 0 => Some((c - 1, self.chunks[c - 1].len() - 1)),
            _ => None,
        };
        (prev, self.first_from(slot))
    }
    /// Removes the entry at `slot` and repairs `letter_changes`, the neighbour LCP and the histogram.
    /// Key fact: the list is sorted, so after removing e the LCP of its old neighbours is
    /// min(lcp(prev,e), lcp(e,next)); no string comparison is needed.
    fn remove(&mut self, slot: Slot) {
        let e = self.at(slot);
        let prev = self.neighbours_of_gap(slot).0;
        let (p, n) = (prev.map(|s| self.at(s)), self.first_from((slot.0, slot.1 + 1)).map(|s| self.at(s)));
        if let Some(p) = p {
            self.letter_changes -= (p.letter != e.letter) as usize; // pair (prev, e) disappears
            self.histogram_remove(p.lcp_next);
        }
        if let Some(n) = n {
            self.letter_changes -= (e.letter != n.letter) as usize; // pair (e, next) disappears
            self.histogram_remove(e.lcp_next);
        }
        let mut joined = 0;
        if let (Some(p), Some(n)) = (p, n) {
            self.letter_changes += (p.letter != n.letter) as usize; // pair (prev, next) appears
            joined = p.lcp_next.min(e.lcp_next);
            self.histogram_add(joined);
        }
        if let Some(p) = prev {
            self.chunks[p.0][p.1].lcp_next = joined;
        }
        self.chunks[slot.0].remove(slot.1);
        if self.chunks[slot.0].is_empty() {
            self.remove_chunk(slot.0);
        }
        self.len -= 1;
    }
    /// Inserts a new entry in front of `slot`. The caller passes the LCP with the previous and the
    /// next entry (they are measured by the caller against the two neighbours).
    fn insert(&mut self, slot: Slot, id: u32, letter: u8, lcp_prev: u32, lcp_next: u32) {
        if self.chunks.is_empty() {
            self.insert_chunk(0, vec![]);
        }
        let (prev, next) = self.neighbours_of_gap(slot);
        let (p, n) = (prev.map(|s| self.at(s)), next.map(|s| self.at(s)));
        let mut entry = Entry { id, letter, lcp_next: 0 };
        if let (Some(ps), Some(p)) = (prev, p) {
            self.letter_changes += (p.letter != letter) as usize;
            self.histogram_add(lcp_prev);
            self.chunks[ps.0][ps.1].lcp_next = lcp_prev;
        }
        if let Some(n) = n {
            self.letter_changes += (letter != n.letter) as usize;
            self.histogram_add(lcp_next);
            entry.lcp_next = lcp_next;
        }
        if let (Some(p), Some(n)) = (p, n) {
            self.letter_changes -= (p.letter != n.letter) as usize; // the old pair (prev, next) is split
            self.histogram_remove(p.lcp_next);
        }
        self.chunks[slot.0].insert(slot.1, entry);
        self.chunk_of[id as usize] = self.chunk_uid[slot.0];
        self.len += 1;
        if self.chunks[slot.0].len() > 2 * CHUNK {
            let second_half = self.chunks[slot.0].split_off(CHUNK); // chunk too big: split in two halves
            self.insert_chunk(slot.0 + 1, second_half);
        }
    }
}

// ------------------------------------------ the sweep ------------------------------------------

/// `table[p]` = length of the longest Lyndon word that starts at `text[p]` and ends inside `text`
/// (the "Lyndon table" of `text`). Computed from the ranks of the rotations: for a Lyndon word the
/// rotation order equals the suffix order, and the longest Lyndon word starting at p ends where
/// the next smaller suffix starts. If your library already builds this table, use it instead.
fn lyndon_table(sorted_ids: impl Iterator<Item = usize>, text_len: usize) -> Vec<usize> {
    let mut rank = vec![0usize; text_len];
    for (r, id) in sorted_ids.enumerate() {
        rank[id] = r;
    }
    let mut table = vec![0usize; text_len];
    let mut stack: Vec<usize> = Vec::new(); // positions with increasing rank, from the right
    for p in (0..text_len).rev() {
        while stack.last().is_some_and(|&q| rank[q] > rank[p]) {
            stack.pop();
        }
        table[p] = stack.last().copied().unwrap_or(text_len) - p; // next smaller suffix - p
        stack.push(p);
    }
    table
}

/// Splits the ranges of `a` and of `b` (each ascending and disjoint) into maximal pieces on which
/// "is in a" and "is in b" do not change, and calls `f(lo, hi, in_a, in_b)` for every piece that is
/// in at least one of them. Costs O(|a| + |b|): no position is visited unless `f` does it.
fn for_each_piece(a: &[Factor], b: &[Factor], mut f: impl FnMut(usize, usize, bool, bool)) {
    let (mut i, mut j, mut cur) = (0, 0, 0);
    while i < a.len() || j < b.len() {
        let (a_lo, a_hi) = a.get(i).map_or((usize::MAX, usize::MAX), |r| (r.start, r.end()));
        let (b_lo, b_hi) = b.get(j).map_or((usize::MAX, usize::MAX), |r| (r.start, r.end()));
        cur = cur.max(a_lo.min(b_lo)); // skip the gap before the next range
        let (in_a, in_b) = (a_lo <= cur, b_lo <= cur);
        let hi = (if in_a { a_hi } else { a_lo }).min(if in_b { b_hi } else { b_lo });
        f(cur, hi, in_a, in_b);
        cur = hi;
        i += (in_a && a_hi == hi) as usize;
        j += (in_b && b_hi == hi) as usize;
    }
}

/// State of the sweep over all rotations of `text`: where we are, the factorization of the current
/// conjugate, and the sorted entries.
#[derive(Clone)]
struct RotationSweep {
    text_len: usize,
    window_start: usize, // the current conjugate is text[window_start..] + text[..window_start]
    doubled: Vec<u8>,       // text·text
    lyndon_len: Vec<usize>, // Lyndon table of `text` (see `lyndon_table`)
    prefix_hash: Vec<u64>, // prefix_hash[k] = polynomial hash of the first k chars of text·text
    base_power: Vec<u64>,  // base_power[k] = HASH_BASE^k mod HASH_MOD
    factors: VecDeque<Factor>, // Lyndon factors of the current conjugate, left to right, tiling [window_start, window_start+text_len)
    sorted: SortedEntries,     // the entries, in omega-sorted order
    pending: Vec<u32>, // ids that must be inserted but have been postponed (may hold stale ids, see is_pending)
    is_pending: Vec<bool>, // is_pending[id]: the id is really pending
    pending_count: usize, // number of ids with is_pending == true
    hint_offset: Vec<u32>, // hint_offset[id] = index inside its chunk that a deleted entry had (u32::MAX: no hint)
    stats: ChainStats,
}

impl RotationSweep {
    /// Builds the starting state: window_start = 0, `text` is Lyndon, so it is ONE factor and its
    /// entries are its text_len rotations. `circular_sa` = starting indices of the rotations of `text` in sorted
    /// order (the circular suffix array); if given it replaces the initial comparison sort.
    fn new(text: &[u8], circular_sa: Option<&[usize]>) -> Self {
        let text_len = text.len();
        let doubled = [text, text].concat();
        let mut prefix_hash = vec![0u64; 2 * text_len + 1];
        for k in 0..2 * text_len {
            prefix_hash[k + 1] = (mul_mod(prefix_hash[k], HASH_BASE) + doubled[k] as u64 + 1) % HASH_MOD;
        }
        let mut base_power = vec![1u64; text_len + 1];
        for k in 0..text_len {
            base_power[k + 1] = mul_mod(base_power[k], HASH_BASE);
        }
        let mut stats = ChainStats::default();
        let rotation = |p: usize| Conjugate { start: 0, len: text_len, off: p };
        let sa: Vec<usize> = match circular_sa {
            Some(sa) => sa.to_vec(),
            None => {
                // no suffix array given: sort the text_len rotations by comparison
                let mut v: Vec<usize> = (0..text_len).collect();
                v.sort_by(|&a, &b| compare_conjugates(&doubled, &rotation(a), &rotation(b), &mut stats).1);
                v
            }
        };
        let lcps: Vec<u32> = sa
            .windows(2)
            .map(|w| compare_conjugates(&doubled, &rotation(w[0]), &rotation(w[1]), &mut stats).0 as u32)
            .collect();
        // the BBWT letter of the rotation starting at p is the character before it
        let entries: Vec<(u32, u8)> = sa.iter().map(|&p| (p as u32, doubled[p + text_len - 1])).collect();
        RotationSweep {
            text_len,
            window_start: 0,
            lyndon_len: lyndon_table(sa.iter().copied(), text_len),
            prefix_hash,
            base_power,
            factors: VecDeque::from(vec![Factor { start: 0, len: text_len }]),
            sorted: SortedEntries::from_sorted(text_len, &entries, &lcps),
            pending: vec![],
            is_pending: vec![false; text_len],
            pending_count: 0,
            hint_offset: vec![u32::MAX; text_len],
            doubled,
            stats,
        }
    }

    /// Hash of the text of one factor occurrence (difference of two prefix hashes).
    fn factor_hash(&self, f: Factor) -> u64 {
        (self.prefix_hash[f.end()] + HASH_MOD - mul_mod(self.prefix_hash[f.start], self.base_power[f.len])) % HASH_MOD
    }
    /// Are the two factor occurrences the same string? (same length and same hash; no char check)
    fn same_factor(&self, f: Factor, g: Factor) -> bool {
        f.len == g.len && self.factor_hash(f) == self.factor_hash(g)
    }
    fn runs(&self) -> usize {
        self.sorted.runs()
    }
    /// The BBWT letter of a conjugate: the character just before it inside its factor
    /// (cyclically: for offset 0 it is the last character of the factor).
    fn predecessor_letter(&self, c: &Conjugate) -> u8 {
        self.doubled[if c.off == 0 { c.start + c.len - 1 } else { c.start + c.off - 1 }]
    }
    fn mark_pending(&mut self, id: usize) {
        if !self.is_pending[id] {
            self.is_pending[id] = true;
            self.pending_count += 1;
            self.pending.push(id as u32);
        }
    }
    fn unmark_pending(&mut self, id: usize) {
        if self.is_pending[id] {
            self.is_pending[id] = false;
            self.pending_count -= 1;
            self.hint_offset[id] = u32::MAX;
        }
    }

    /// LCP of `conj` with the entry at `slot` (0 if there is none).
    fn lcp_with(&mut self, slot: Option<Slot>, conj: &Conjugate) -> u32 {
        let Some(slot) = slot else { return 0 };
        let other = conjugate_of(self.text_len, self.window_start, self.sorted.at(slot).id, &self.factors);
        compare_conjugates(&self.doubled, &other, conj, &mut self.stats).0 as u32
    }

    /// Inserts entry `id` into the sorted list (current factors): search for its place (starting from
    /// the index it had before, if it was deleted this way), then insert with the LCPs measured
    /// against the two neighbours.
    fn insert_entry(&mut self, id: u32) {
        let (text_len, window_start) = (self.text_len, self.window_start);
        let conj = conjugate_of(text_len, window_start, id, &self.factors);
        let hint = std::mem::replace(&mut self.hint_offset[id as usize], u32::MAX);
        let slot = {
            let Self { doubled, factors, sorted, stats, .. } = self;
            let compare_to_target = |other: u32| {
                compare_conjugates(doubled, &conjugate_of(text_len, window_start, other, factors), &conj, stats).1
            };
            if hint == u32::MAX {
                sorted.seek(compare_to_target)
            } else {
                sorted.seek_near(id, hint as usize, compare_to_target)
            }
        };
        let (prev, next) = self.sorted.neighbours_of_gap(slot);
        let (lcp_prev, lcp_next) = (self.lcp_with(prev, &conj), self.lcp_with(next, &conj));
        let letter = self.predecessor_letter(&conj);
        self.sorted.insert(slot, id, letter, lcp_prev, lcp_next);
    }

    /// Inserts every pending id (skipping ids that were cancelled in the meantime).
    fn flush_pending(&mut self) {
        if self.pending_count > 0 {
            self.stats.flushes += 1;
            for k in 0..self.pending.len() {
                let id = self.pending[k];
                if std::mem::replace(&mut self.is_pending[id as usize], false) {
                    self.insert_entry(id);
                    self.stats.phys_ins += 1;
                }
            }
            self.pending_count = 0;
        }
        self.pending.clear();
    }

    /// Takes entry `id` out of the sorted list. If `requeue`, the entry is inserted again later with
    /// its new word, and the index it had is remembered as the starting point of that search.
    /// An id that is already pending is not in the list: nothing to remove (this also makes a
    /// second visit of the same risky id harmless).
    fn delete_entry(&mut self, id: usize, requeue: bool) {
        if self.is_pending[id] {
            if !requeue {
                self.unmark_pending(id);
            }
            return;
        }
        let slot = self.sorted.find(id as u32);
        self.sorted.remove(slot);
        self.stats.phys_del += 1;
        if requeue {
            self.stats.risky += 1;
            self.mark_pending(id);
            self.hint_offset[id] = slot.1 as u32;
        }
    }

    /// One left rotation: window_start -> window_start + 1. Updates the factors and the sorted list.
    /// `best` = fewest runs found so far (used to decide whether the re-insertions can be skipped).
    fn rotate_left_once(&mut self, prune: bool, best: usize) {
        let (text_len, old_start) = (self.text_len, self.window_start);
        self.stats.steps += 1;
        let max_lcp = self.sorted.max_lcp; // read before any deletion changes it
        self.stats.lmax_sum += max_lcp as u64;

        // ---- A. new factorization: split the first factor, append the moved char, merge at the back ----
        let old_count = self.factors.len();
        let first = self.factors.pop_front().unwrap();
        // The first block survives iff it has a second copy: a merged factor starts with a letter
        // smaller than the moved char, the copy starts with that char (I1), so it is never merged.
        let first_block_survives = old_count > 1 && self.same_factor(first, self.factors[0]);
        // pieces = the Lyndon factors of first[1..], by walking the Lyndon table (I2)
        let mut pieces = Vec::new();
        let mut p = first.start + 1;
        while p < first.end() {
            let len = self.lyndon_len[p];
            pieces.push(Factor { start: p, len });
            p += len;
        }
        for &piece in pieces.iter().rev() {
            self.factors.push_front(piece);
        }
        // the moved char becomes a new one-char factor at the back, then merges with the factors before it
        let mut back = Factor { start: old_start + text_len, len: 1 };
        self.factors.push_back(back);
        let mut merged_owners = Vec::new(); // the last copy of each merged block, from the back
        let mut right_merged: Option<Factor> = None; // the factor merged just before, i.e. to the right
        let mut old_left = old_count - 1; // old factors between the pieces and `back`
        while old_left > 0 {
            let before = self.factors[self.factors.len() - 2];
            if !factor_is_smaller(&self.doubled, before, back) {
                break;
            }
            self.factors.pop_back();
            if right_merged.map_or(true, |r| !self.same_factor(before, r)) {
                merged_owners.push(before);
            }
            right_merged = Some(before);
            old_left -= 1;
            back = Factor { start: before.start, len: before.len + back.len };
            *self.factors.back_mut().unwrap() = back;
        }
        self.window_start = old_start + 1;

        // ---- B. which owners die, which are born (both lists ascending, disjoint) ----
        // Dead owners: the first factor without a second copy, and every merged block (by its last
        // copy). The first position of the window leaves at the front and re-enters at the back,
        // so the first factor's positions are [old_start+1, first.end) and old_start+text_len.
        let dead_first = !first_block_survives;
        let mut dead_owners = Vec::new();
        if dead_first && first.len > 1 {
            dead_owners.push(Factor { start: old_start + 1, len: first.len - 1 });
        }
        dead_owners.extend(merged_owners.iter().rev());
        if dead_first {
            dead_owners.push(Factor { start: old_start + text_len, len: 1 });
        }
        // If `back` equals the old factor before it, no new factor is born: that block just grows,
        // and its entries (already in the list) are renamed to the new last copy, see E2.
        let grown_block = (old_left > 0)
            .then(|| self.factors[self.factors.len() - 2])
            .filter(|&before_back| self.same_factor(before_back, back));
        self.stats.growth += grown_block.is_some() as u64;
        // Born owners: the last piece of each block of equal pieces, and `back` unless it only extends a block.
        let mut born_owners = Vec::new();
        for k in 0..pieces.len() {
            if k + 1 == pieces.len() || !self.same_factor(pieces[k], pieces[k + 1]) {
                born_owners.push(pieces[k]);
            }
        }
        if grown_block.is_none() {
            born_owners.push(back);
        }

        // ---- C, D, E. one pass over the positions of the owners ----
        for_each_piece(&dead_owners, &born_owners, |lo, hi, in_dead, in_born| match (in_dead, in_born) {
            // C. only dead: the entry exists only before the rotation; delete it for good
            (true, false) => {
                (lo..hi).for_each(|q| self.delete_entry(wrap(q, text_len), false));
                self.stats.dead += (hi - lo) as u64;
            }
            // D. dead and born ("kept"): the entry may change place. Its word changes within max_lcp
            // characters of the end of its old or new factor, i.e. of the end of this piece, and its
            // letter changes at the first position of a factor, i.e. of this piece. Such entries are
            // deleted now and inserted again later (in the new order).
            (true, true) => {
                let tail = lo.max(hi.saturating_sub(max_lcp));
                (tail > lo).then_some(lo).into_iter().chain(tail..hi).for_each(|q| self.delete_entry(wrap(q, text_len), true));
            }
            // E. only born: the entry exists only after the rotation; queue it for insertion
            _ => {
                (lo..hi).for_each(|q| self.mark_pending(wrap(q, text_len)));
                self.stats.ins += (hi - lo) as u64;
            }
        });

        // ---- E2. a grown block: its entries move to the new last copy, ids shift by |back| ----
        // The order of the entries does not change (same conjugates of the same factor). The new
        // ids are free: the positions of `back` were deleted above.
        if let Some(old_last) = grown_block {
            for k in 0..back.len {
                let (old_id, new_id) = (wrap(old_last.start + k, text_len), wrap(back.start + k, text_len));
                if self.is_pending[old_id] {
                    self.unmark_pending(old_id);
                    self.mark_pending(new_id);
                } else {
                    let slot = self.sorted.find(old_id as u32);
                    self.sorted.rename(slot, new_id as u32);
                }
            }
        }

        // ---- F. insert, unless this conjugate provably cannot beat `best` ----
        if prune && self.pending_count > 0 && self.runs() >= best {
            self.stats.deferred_steps += 1;
        } else {
            self.flush_pending();
        }
    }
}

// ------------------------------- oracle & invariants (verification) -------------------------------

fn lyndon_factors_plain(x: &[u8]) -> Vec<(usize, usize)> {
    let (mut out, mut i) = (Vec::new(), 0);
    while i < x.len() {
        let (mut j, mut k) = (i + 1, i);
        while j < x.len() && x[k] <= x[j] {
            k = if x[k] < x[j] { i } else { k + 1 };
            j += 1;
        }
        let l = j - k;
        while i <= k {
            out.push((i, l));
            i += l;
        }
    }
    out
}

/// Independent slow reference: the BBWT letters of x built from scratch (sort all conjugates of the
/// distinct Lyndon factors with explicit infinite-word comparison).
pub fn oracle_letters(x: &[u8]) -> Vec<u8> {
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let mut items: Vec<(Vec<u8>, u8)> = Vec::new();
    for (start, len) in lyndon_factors_plain(x) {
        let f = x[start..start + len].to_vec();
        if !seen.insert(f.clone()) {
            continue;
        }
        for k in 0..len {
            let mut c = f[k..].to_vec();
            c.extend_from_slice(&f[..k]);
            items.push((c, f[(k + len - 1) % len]));
        }
    }
    items.sort_by(|a, b| {
        let l = a.0.len() + b.0.len();
        a.0.iter().cycle().take(l).cmp(b.0.iter().cycle().take(l))
    });
    items.into_iter().map(|e| e.1).collect()
}

pub fn runs_of_letters(v: &[u8]) -> usize {
    (!v.is_empty()) as usize + v.windows(2).filter(|w| w[0] != w[1]).count()
}

impl RotationSweep {
    fn letters(&self) -> Vec<u8> {
        self.sorted.chunks.iter().flatten().map(|e| e.letter).collect()
    }
    /// Checks that the sorted list is really sorted, that every stored letter/LCP/counter is right.
    fn check_invariants(&mut self) -> Result<(), String> {
        let flat: Vec<Entry> = self.sorted.chunks.iter().flatten().copied().collect();
        if flat.len() != self.sorted.len {
            return Err("len".into());
        }
        if self.sorted.chunks.iter().any(|c| c.is_empty() || c.len() > 2 * CHUNK) {
            return Err("chunk size".into());
        }
        for (i, chunk) in self.sorted.chunks.iter().enumerate() {
            let uid = self.sorted.chunk_uid[i];
            if self.sorted.chunk_index[uid as usize] != i || chunk.iter().any(|e| self.sorted.chunk_of[e.id as usize] != uid) {
                return Err("chunk_index / chunk_of".into());
            }
        }
        let mut scratch_stats = ChainStats::default();
        let mut max_lcp = 0;
        for j in 0..flat.len() {
            let conj = conjugate_of(self.text_len, self.window_start, flat[j].id, &self.factors);
            if self.predecessor_letter(&conj) != flat[j].letter {
                return Err(format!("letter at {}", j));
            }
            if j + 1 < flat.len() {
                let next = conjugate_of(self.text_len, self.window_start, flat[j + 1].id, &self.factors);
                let (lcp, order) = compare_conjugates(&self.doubled, &conj, &next, &mut scratch_stats);
                if order != Ordering::Less {
                    return Err(format!("order at {}", j));
                }
                if flat[j].lcp_next as usize != lcp {
                    return Err(format!("lcp at {}", j));
                }
                max_lcp = max_lcp.max(lcp);
            } else if flat[j].lcp_next != 0 {
                return Err("last lcp".into());
            }
        }
        if max_lcp != self.sorted.max_lcp {
            return Err(format!("max_lcp {} vs {}", self.sorted.max_lcp, max_lcp));
        }
        if flat.windows(2).filter(|w| w[0].letter != w[1].letter).count() != self.sorted.letter_changes {
            return Err("letter_changes".into());
        }
        Ok(())
    }
}

/// Sweeps all conjugates text[i..]+text[..i] of the Lyndon word `text`. Returns (index of the first
/// conjugate with the fewest BBWT runs, that run count, statistics). `verify` cross-checks every step
/// against the slow oracle.
pub fn chain_min_runs(text: &[u8], prune: bool, verify: bool) -> (usize, usize, ChainStats) {
    chain_min_runs_sa(text, None, prune, verify)
}

/// Same, with the circular suffix array of `text` (rotation starts in sorted order) supplied by the caller.
pub fn chain_min_runs_sa(
    text: &[u8],
    circular_sa: Option<&[usize]>,
    prune: bool,
    verify: bool,
) -> (usize, usize, ChainStats) {
    let text_len = text.len();
    if text_len < 2 {
        return (0, text_len, ChainStats::default());
    }
    let mut sweep = RotationSweep::new(text, circular_sa);
    let (mut best, mut best_shift) = (sweep.runs(), 0);
    if verify {
        assert_eq!(sweep.letters(), oracle_letters(text), "initial state");
        let table: Vec<usize> = (0..text_len).map(|p| lyndon_factors_plain(&text[p..])[0].1).collect();
        assert_eq!(sweep.lyndon_len, table, "Lyndon table");
    }
    for _ in 0..text_len - 1 {
        let best_before = best;
        sweep.rotate_left_once(prune, best);
        // the run count is exact only when nothing is pending
        if sweep.pending_count == 0 && sweep.runs() < best {
            best = sweep.runs();
            best_shift = sweep.window_start;
        }
        if verify {
            let x: Vec<u8> = (0..text_len).map(|k| text[(sweep.window_start + k) % text_len]).collect();
            let expected = oracle_letters(&x);
            let true_runs = runs_of_letters(&expected);
            let mut flushed = sweep.clone();
            flushed.flush_pending();
            flushed
                .check_invariants()
                .unwrap_or_else(|e| panic!("invariant after step to {}: {}", sweep.window_start, e));
            assert_eq!(
                flushed.letters(),
                expected,
                "letters differ at rotation {} (pending {})",
                sweep.window_start,
                sweep.pending_count
            );
            if sweep.pending_count > 0 {
                assert!(sweep.runs() <= true_runs, "partial list has more runs than the true count");
                assert!(
                    true_runs >= best_before,
                    "deferred a conjugate that beats the best: {} < {}",
                    true_runs,
                    best_before
                );
            } else {
                assert_eq!(sweep.runs(), true_runs);
            }
        }
    }
    sweep.flush_pending();
    (best_shift, best, sweep.stats.clone())
}

// ------------------------------------- public entry points -------------------------------------

/// The least rotation of the primitive root of `input`: the Lyndon word whose conjugates are swept.
fn lyndon_root(input: &[u8]) -> Vec<u8> {
    let (start, period) = shiloach(input);
    (0..period).map(|k| input[(start + k) % input.len()]).collect()
}

/// Index (among the conjugates of the least rotation of the primitive root) of the first conjugate
/// with the fewest BBWT runs. `prune` enables the certified skipping.
pub fn optimal_bbwt_chain(input: &[u8], prune: bool) -> usize {
    let text = lyndon_root(input);
    let circular_sa = circular_suffix_array(&text); // initial order of the text_len rotations: the only suffix-array computation of the run
    chain_min_runs_sa(&text, Some(&circular_sa), prune, false).0
}

/// Same sweep, cross-checked against the oracle at every step when n <= 400, with statistics.
pub fn online_bbwt_check(input: &[u8]) -> usize {
    let text = lyndon_root(input);
    let verify = text.len() <= 400;
    let circular_sa = circular_suffix_array(&text);
    let (shift, runs, st) = chain_min_runs_sa(&text, Some(&circular_sa), true, verify);
    let steps = st.steps.max(1) as f64;
    println!("n={} verified={} best shift={} runs={} | per step: risky {:.1} deleted {:.1} reinserted {:.1} cmps {:.0} | deferred {} of {} steps, block growths {}",
        text.len(), verify, shift, runs, st.risky as f64 / steps, st.phys_del as f64 / steps, st.phys_ins as f64 / steps,
        st.cmp_calls as f64 / steps, st.deferred_steps, st.steps, st.growth);
    shift
}
