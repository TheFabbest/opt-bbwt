use crate::lyndon::lyndon_table;
use crate::lyndon::shiloach;
use crate::suffix_array::{circular_suffix_array, sa_to_rank};
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
//             Positions are "absolute": window_start..window_start+text_len; `wrap_index` maps them to 0..text_len.
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
//   Sorted entries live in chunks of ~512 (insert/delete cost O(chunk)). Each entry stores its LCP
//   with the next entry, so deleting an entry gives min(lcp_prev, lcp_self) for the joined pair,
//   and a histogram of all LCPs gives max_lcp in O(1).
// =============================================================================================

const HASH_MOD: u64 = (1u64 << 61) - 1;
const HASH_BASE: u64 = 1_000_003_117_927_911;
const CHUNK: usize = 512;

/// (a * b) mod (2^61 - 1): 128-bit product, then fold the high bits onto the low ones.
fn mul_mod(a: u64, b: u64) -> u64 {
    let product = (a as u128) * (b as u128);
    let mut r = (product & HASH_MOD as u128) as u64 + (product >> 61) as u64;
    if r >= HASH_MOD {
        r -= HASH_MOD;
    }
    r
}

/// One occurrence of a Lyndon factor in the window: text[start..start+len] (absolute positions).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Factor {
    start: usize,
    len: usize,
}
impl Factor {
    fn end(&self) -> usize {
        self.start + self.len
    }
}

/// Returns x mod text_len, for 0 <= x < 3*text_len, using comparisons instead of a division.
/// It maps an absolute position (which may run past the end of the text because the window wraps
/// around) back to an index of the text. It is called in the innermost comparison loop, where a
/// real `%` was measurably slower.
#[inline(always)]
fn wrap_index(x: usize, text_len: usize) -> usize {
    if x < text_len {
        x
    } else if x < 2 * text_len {
        x - text_len
    } else {
        x - 2 * text_len
    }
}

/// A conjugate of a factor: the factor's text rotated left by `offset`. It is infinite-periodic
/// when compared (see `compare_conjugates`).
#[derive(Clone, Copy)]
struct Conjugate {
    factor_start: usize,
    factor_len: usize,
    offset: usize,
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

/// Compares two conjugates as infinite periodic words (each one repeats its own factor forever).
/// Returns (length of their longest common prefix, Less/Equal/Greater).
/// How it works: walk both words in parallel, each inside its own factor, wrapping to the factor's
/// start when it ends. Two periodic words with periods la, lb that agree on la+lb characters are
/// equal (Fine–Wilf), so `la + lb` characters are enough to decide.
fn compare_conjugates(
    text: &[u8],
    text_len: usize,
    a: &Conjugate,
    b: &Conjugate,
    stats: &mut ChainStats,
) -> (usize, Ordering) {
    stats.cmp_calls += 1;
    let max_len = a.factor_len + b.factor_len;
    let (mut off_a, mut off_b) = (a.offset, b.offset); // offset inside each factor
    let (base_a, base_b) = (
        wrap_index(a.factor_start, text_len),
        wrap_index(b.factor_start, text_len),
    );
    for k in 0..max_len {
        let mut pos_a = base_a + off_a;
        if pos_a >= text_len {
            pos_a -= text_len;
        } // text index of the k-th char of a
        let mut pos_b = base_b + off_b;
        if pos_b >= text_len {
            pos_b -= text_len;
        }
        let (char_a, char_b) = (text[pos_a], text[pos_b]);
        if char_a != char_b {
            stats.cmp_chars += (k + 1) as u64;
            return (k, char_a.cmp(&char_b));
        }
        off_a += 1;
        if off_a == a.factor_len {
            off_a = 0;
        } // wrap inside the factor
        off_b += 1;
        if off_b == b.factor_len {
            off_b = 0;
        }
    }
    stats.cmp_chars += max_len as u64;
    (max_len, Ordering::Equal)
}

/// The conjugate that entry `id` stands for, in the window starting at `window_start`
/// whose factors are `factors`. (id -> absolute position -> the factor containing it.)
fn conjugate_of_entry(
    text_len: usize,
    id: u32,
    factors: &VecDeque<Factor>,
    window_start: usize,
) -> Conjugate {
    conjugate_at(
        window_start
            + wrap_index(
                id as usize + text_len - wrap_index(window_start, text_len),
                text_len,
            ),
        factors,
    )
}

/// The conjugate that starts at absolute position `abs`: finds the factor containing `abs`
/// by binary search on the sorted factor starts.
#[inline(always)]
fn conjugate_at(abs: usize, factors: &VecDeque<Factor>) -> Conjugate {
    let f = factors[factors.partition_point(|f| f.start <= abs) - 1];
    Conjugate {
        factor_start: f.start,
        factor_len: f.len,
        offset: abs - f.start,
    }
}

/// Same as `conjugate_of_entry` but for the window BEFORE this rotation (it started at
/// `old_window_start`). `factors` is already the NEW factor list, so the old factors are rebuilt from
/// `removed` (the old factors that no longer exist): if the position lies inside a removed factor,
/// that factor is the old one; otherwise the position was not touched and the new list is right.
fn conjugate_of_entry_before_rotation(
    text_len: usize,
    old_window_start: usize,
    id: u32,
    removed: &[Factor],
    factors: &VecDeque<Factor>,
) -> Conjugate {
    let abs = old_window_start
        + wrap_index(
            id as usize + text_len - wrap_index(old_window_start, text_len),
            text_len,
        );
    let k = removed.partition_point(|f| f.start <= abs);
    if k > 0 {
        let f = removed[k - 1];
        if abs < f.start + f.len {
            return Conjugate {
                factor_start: f.start,
                factor_len: f.len,
                offset: abs - f.start,
            };
        }
    }
    conjugate_at(abs, factors) // abs is in [old_window_start+1, old_window_start+text_len): same factor in old and new window
}

/// Compares the two factors f and g as plain strings in lexicographic order (a proper prefix is
/// smaller) and returns f < g. Used for the Lyndon merge: if f < g then f·g is a Lyndon word,
/// so two neighbouring factors f, g with f < g must be merged into one.
/// How it works: compare character by character up to the shorter length; the first difference
/// decides; if there is none, the shorter one is smaller.
fn factor_is_smaller(text: &[u8], text_len: usize, f: Factor, g: Factor) -> bool {
    for k in 0..f.len.min(g.len) {
        let (cf, cg) = (
            text[wrap_index(f.start + k, text_len)],
            text[wrap_index(g.start + k, text_len)],
        );
        if cf != cg {
            return cf < cg;
        }
    }
    f.len < g.len
}

// ------------------------------------ sorted entry list ------------------------------------

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
    len: usize,
    letter_changes: usize, // number of adjacent pairs with different letters; runs = this + 1
    lcp_histogram: Vec<u32>, // lcp_histogram[v] = how many adjacent pairs have LCP exactly v
    max_lcp: usize,        // largest v with lcp_histogram[v] > 0
}

impl SortedEntries {
    fn new(text_len: usize) -> Self {
        SortedEntries {
            chunks: vec![],
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
        list.chunks = all.chunks(CHUNK).map(|c| c.to_vec()).collect();
        list
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
    /// The slots just before and just after the gap in front of `slot` (None at the list ends).
    /// `slot.1` may equal the chunk length, meaning "after the last element of the chunk".
    fn neighbours_of_gap(&self, slot: Slot) -> (Option<Slot>, Option<Slot>) {
        let prev = if slot.1 > 0 {
            Some((slot.0, slot.1 - 1))
        } else if slot.0 > 0 {
            Some((slot.0 - 1, self.chunks[slot.0 - 1].len() - 1))
        } else {
            None
        };
        let next = if slot.1 < self.chunks[slot.0].len() {
            Some(slot)
        } else if slot.0 + 1 < self.chunks.len() {
            Some((slot.0 + 1, 0))
        } else {
            None
        };
        (prev, next)
    }
    /// Removes the entry at `slot` and repairs `letter_changes`, the neighbour LCP and the histogram.
    /// Key fact: the list is sorted, so after removing e the LCP of its old neighbours is
    /// min(lcp(prev,e), lcp(e,next)); no string comparison is needed.
    fn remove(&mut self, slot: Slot) {
        let removed = self.at(slot);
        let prev = self.neighbours_of_gap(slot).0;
        let next = if slot.1 + 1 < self.chunks[slot.0].len() {
            Some((slot.0, slot.1 + 1))
        } else if slot.0 + 1 < self.chunks.len() {
            Some((slot.0 + 1, 0))
        } else {
            None
        };
        match (prev, next) {
            (Some(prev_slot), Some(next_slot)) => {
                let (prev_entry, next_entry) = (self.at(prev_slot), self.at(next_slot));
                if prev_entry.letter != removed.letter {
                    self.letter_changes -= 1;
                } // pair (prev, removed) disappears
                if removed.letter != next_entry.letter {
                    self.letter_changes -= 1;
                } // pair (removed, next) disappears
                if prev_entry.letter != next_entry.letter {
                    self.letter_changes += 1;
                } // pair (prev, next) appears
                let joined_lcp = prev_entry.lcp_next.min(removed.lcp_next);
                self.histogram_remove(prev_entry.lcp_next);
                self.histogram_remove(removed.lcp_next);
                self.histogram_add(joined_lcp);
                self.chunks[prev_slot.0][prev_slot.1].lcp_next = joined_lcp;
            }
            (Some(prev_slot), None) => {
                // removing the last entry
                let prev_entry = self.at(prev_slot);
                if prev_entry.letter != removed.letter {
                    self.letter_changes -= 1;
                }
                self.histogram_remove(prev_entry.lcp_next);
                self.chunks[prev_slot.0][prev_slot.1].lcp_next = 0;
            }
            (None, Some(next_slot)) => {
                // removing the first entry
                if removed.letter != self.at(next_slot).letter {
                    self.letter_changes -= 1;
                }
                self.histogram_remove(removed.lcp_next);
            }
            (None, None) => {}
        }
        self.chunks[slot.0].remove(slot.1);
        if self.chunks[slot.0].is_empty() {
            self.chunks.remove(slot.0);
        }
        self.len -= 1;
    }
    /// Inserts a new entry in front of `slot`. The caller passes the LCP with the previous and the
    /// next entry (they were computed while comparing during the binary search).
    fn insert(&mut self, slot: Slot, id: u32, letter: u8, lcp_with_prev: u32, lcp_with_next: u32) {
        if self.chunks.is_empty() {
            self.chunks.push(vec![Entry {
                id,
                letter,
                lcp_next: 0,
            }]);
            self.len = 1;
            return;
        }
        let (prev, next) = self.neighbours_of_gap(slot);
        let mut new_entry = Entry {
            id,
            letter,
            lcp_next: 0,
        };
        match (prev, next) {
            (Some(prev_slot), Some(next_slot)) => {
                let (prev_entry, next_entry) = (self.at(prev_slot), self.at(next_slot));
                if prev_entry.letter != letter {
                    self.letter_changes += 1;
                }
                if letter != next_entry.letter {
                    self.letter_changes += 1;
                }
                if prev_entry.letter != next_entry.letter {
                    self.letter_changes -= 1;
                } // the old pair (prev, next) is split
                self.histogram_remove(prev_entry.lcp_next);
                self.histogram_add(lcp_with_prev);
                self.histogram_add(lcp_with_next);
                self.chunks[prev_slot.0][prev_slot.1].lcp_next = lcp_with_prev;
                new_entry.lcp_next = lcp_with_next;
            }
            (Some(prev_slot), None) => {
                // inserting at the end
                if self.at(prev_slot).letter != letter {
                    self.letter_changes += 1;
                }
                self.histogram_add(lcp_with_prev);
                self.chunks[prev_slot.0][prev_slot.1].lcp_next = lcp_with_prev;
            }
            (None, Some(next_slot)) => {
                // inserting at the front
                if letter != self.at(next_slot).letter {
                    self.letter_changes += 1;
                }
                self.histogram_add(lcp_with_next);
                new_entry.lcp_next = lcp_with_next;
            }
            (None, None) => unreachable!(),
        }
        self.chunks[slot.0].insert(slot.1, new_entry);
        self.len += 1;
        if self.chunks[slot.0].len() > 2 * CHUNK {
            // chunk too big: split in two halves
            let second_half = self.chunks[slot.0].split_off(CHUNK);
            self.chunks.insert(slot.0 + 1, second_half);
        }
    }
}

// ------------------------------------------ the sweep ------------------------------------------

/// Calls `f(q)` for every position q in [lo, hi) that is covered by none of `owners`
/// (ascending, disjoint). Costs O(#owners overlapping the range + #positions reported), so a
/// long range that is almost entirely covered is never walked position by position.
fn for_each_uncovered(lo: usize, hi: usize, owners: &[Factor], mut f: impl FnMut(usize)) {
    let mut cur = lo;
    for o in &owners[owners.partition_point(|o| o.end() <= lo)..] {
        if o.start >= hi {
            break;
        }
        (cur..o.start).for_each(&mut f);
        cur = cur.max(o.end());
    }
    (cur..hi).for_each(&mut f);
}

/// State of the sweep over all rotations of `text`: where we are, the factorization of the current
/// conjugate, and the sorted entries.
#[derive(Clone)]
struct RotationSweep<'a> {
    text: &'a [u8],
    text_len: usize,
    window_start: usize, // the current conjugate is text[window_start..] + text[..window_start]
    lyndon_len: Vec<usize>, // Lyndon table of `text` (see `lyndon_table`)
    prefix_hash: Vec<u64>, // prefix_hash[k] = polynomial hash of the first k chars of text·text·text
    base_power: Vec<u64>,  // base_power[k] = HASH_BASE^k mod HASH_MOD
    factors: VecDeque<Factor>, // Lyndon factors of the current conjugate, left to right, tiling [window_start, window_start+text_len)
    sorted: SortedEntries,     // the entries, in omega-sorted order
    pending: Vec<u32>, // ids that must be inserted but have been postponed (may hold stale ids, see is_pending)
    is_pending: Vec<bool>, // is_pending[id]: the id is really pending
    pending_count: usize, // number of ids with is_pending == true
    // buffers reused by every rotation (so that a step allocates nothing)
    removed: Vec<Factor>, // old factors that no longer exist, ascending: the first factor, then the merged ones
    pieces: Vec<Factor>,  // factors of the first factor without its first character
    dead_owners: Vec<Factor>, // ascending; see HOW ENTRIES CHANGE
    born_owners: Vec<Factor>, // ascending
    renames: Vec<(usize, usize, Option<Slot>)>, // (old id, new id, slot of the entry) of a grown block
    stats: ChainStats,
}

impl<'a> RotationSweep<'a> {
    /// Builds the starting state: window_start = 0, `text` is Lyndon, so it is ONE factor and its
    /// entries are its text_len rotations. `circular_sa` = starting indices of the rotations of `text` in sorted
    /// order (the circular suffix array); if given it replaces the initial comparison sort.
    fn new(text: &'a [u8], circular_sa: Option<&[usize]>) -> Self {
        let text_len = text.len();
        let total = 3 * text_len + 2;
        let (mut prefix_hash, mut base_power) = (vec![0u64; total + 1], vec![1u64; total + 1]);
        for k in 0..total {
            prefix_hash[k + 1] =
                (mul_mod(prefix_hash[k], HASH_BASE) + text[k % text_len] as u64 + 1) % HASH_MOD;
            base_power[k + 1] = mul_mod(base_power[k], HASH_BASE);
        }
        let mut sweep = RotationSweep {
            text,
            text_len,
            window_start: 0,
            lyndon_len: Vec::new(),
            prefix_hash,
            base_power,
            factors: VecDeque::from(vec![Factor {
                start: 0,
                len: text_len,
            }]),
            sorted: SortedEntries::new(text_len),
            pending: vec![],
            is_pending: vec![false; text_len],
            pending_count: 0,
            removed: vec![],
            pieces: vec![],
            dead_owners: vec![],
            born_owners: vec![],
            renames: vec![],
            stats: ChainStats::default(),
        };
        // (conjugate, id) of the text_len entries, in sorted order
        let sorted_entries: Vec<(Conjugate, u32)> = match circular_sa {
            Some(sa) => sa
                .iter()
                .map(|&p| {
                    (
                        conjugate_of_entry(text_len, p as u32, &sweep.factors, 0),
                        p as u32,
                    )
                })
                .collect(),
            None => {
                // no suffix array given: sort the text_len rotations by comparison
                let mut v: Vec<(Conjugate, u32)> = (0..text_len as u32)
                    .map(|p| (conjugate_of_entry(text_len, p, &sweep.factors, 0), p))
                    .collect();
                let stats = &mut sweep.stats;
                v.sort_by(|a, b| compare_conjugates(text, text_len, &a.0, &b.0, stats).1);
                v
            }
        };
        sweep.lyndon_len = lyndon_table(
            text,
            &sa_to_rank(
                &sorted_entries
                    .iter()
                    .map(|(_, id)| *id as usize)
                    .collect::<Vec<_>>(),
            ),
        );
        let lcps: Vec<u32> = (0..text_len - 1)
            .map(|j| {
                compare_conjugates(
                    text,
                    text_len,
                    &sorted_entries[j].0,
                    &sorted_entries[j + 1].0,
                    &mut sweep.stats,
                )
                .0 as u32
            })
            .collect();
        let list: Vec<(u32, u8)> = sorted_entries
            .iter()
            .map(|(conj, id)| (*id, sweep.predecessor_letter(conj)))
            .collect();
        sweep.sorted = SortedEntries::from_sorted(text_len, &list, &lcps);
        sweep
    }

    /// Hash of the text of one factor occurrence (difference of two prefix hashes).
    fn factor_hash(&self, f: Factor) -> u64 {
        (self.prefix_hash[f.start + f.len] + HASH_MOD
            - mul_mod(self.prefix_hash[f.start], self.base_power[f.len]))
            % HASH_MOD
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
        self.text[wrap_index(
            if c.offset == 0 {
                c.factor_start + c.factor_len - 1
            } else {
                c.factor_start + c.offset - 1
            },
            self.text_len,
        )]
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
        }
    }

    /// Inserts entry `id` into the sorted list (current factors): binary search for its place,
    /// then insert with the LCPs measured against the two neighbours.
    fn insert_entry(&mut self, id: u32) {
        let (text_len, window_start, text) = (self.text_len, self.window_start, self.text);
        let conj = conjugate_of_entry(text_len, id, &self.factors, window_start);
        let slot = {
            let (factors, stats) = (&self.factors, &mut self.stats);
            self.sorted.seek(|other_id| {
                compare_conjugates(
                    text,
                    text_len,
                    &conjugate_of_entry(text_len, other_id, factors, window_start),
                    &conj,
                    stats,
                )
                .1
            })
        };
        let letter = self.predecessor_letter(&conj);
        if self.sorted.chunks.is_empty() {
            self.sorted.insert(slot, id, letter, 0, 0);
            return;
        }
        let (prev, next) = self.sorted.neighbours_of_gap(slot);
        let lcp_with_prev = match prev {
            Some(prev_slot) => {
                let prev_conj = conjugate_of_entry(
                    text_len,
                    self.sorted.at(prev_slot).id,
                    &self.factors,
                    window_start,
                );
                compare_conjugates(text, text_len, &prev_conj, &conj, &mut self.stats).0 as u32
            }
            None => 0,
        };
        let lcp_with_next = match next {
            Some(next_slot) => {
                let next_conj = conjugate_of_entry(
                    text_len,
                    self.sorted.at(next_slot).id,
                    &self.factors,
                    window_start,
                );
                compare_conjugates(text, text_len, &conj, &next_conj, &mut self.stats).0 as u32
            }
            None => 0,
        };
        self.sorted
            .insert(slot, id, letter, lcp_with_prev, lcp_with_next);
    }

    /// Inserts every pending id (skipping ids that were cancelled in the meantime).
    fn flush_pending(&mut self) {
        if self.pending_count == 0 {
            return;
        }
        self.stats.flushes += 1;
        for id in std::mem::take(&mut self.pending) {
            if !self.is_pending[id as usize] {
                continue;
            }
            self.is_pending[id as usize] = false;
            self.insert_entry(id);
            self.stats.phys_ins += 1;
        }
        self.pending_count = 0;
    }

    /// Finds the slot of entry `id` in the sorted list as it was BEFORE this rotation
    /// (the list has not been modified yet; only `factors` has). `removed` = old factors that
    /// vanished. The comparison treats `id` itself as Equal, so the search stops on it and never
    /// compares an entry with itself (that would cost a whole factor length).
    fn locate_entry_before_rotation(
        &mut self,
        old_window_start: usize,
        id: usize,
        removed: &[Factor],
    ) -> Slot {
        let (text_len, text) = (self.text_len, self.text);
        let conj = conjugate_of_entry_before_rotation(
            text_len,
            old_window_start,
            id as u32,
            removed,
            &self.factors,
        );
        let (factors, stats) = (&self.factors, &mut self.stats);
        let (chunk, offset) = self.sorted.seek(|other| {
            if other == id as u32 {
                Ordering::Equal
            } else {
                compare_conjugates(
                    text,
                    text_len,
                    &conjugate_of_entry_before_rotation(
                        text_len,
                        old_window_start,
                        other,
                        removed,
                        factors,
                    ),
                    &conj,
                    stats,
                )
                .1
            }
        });
        let slot = if offset == self.sorted.chunks[chunk].len() {
            (chunk + 1, 0)
        } else {
            (chunk, offset)
        };
        assert!(self.sorted.at(slot).id == id as u32, "entry not found");
        slot
    }

    /// Takes entry `id` out of the sorted list (searching in the list as it was before this
    /// rotation). If `requeue`, the entry is inserted again later with its new word.
    /// An id that is already pending is not in the list: nothing to remove (this also makes a
    /// second visit of the same risky id harmless).
    fn delete_entry(
        &mut self,
        id: usize,
        requeue: bool,
        old_window_start: usize,
        removed: &[Factor],
    ) {
        if self.is_pending[id] {
            if !requeue {
                self.unmark_pending(id);
            }
            return;
        }
        let slot = self.locate_entry_before_rotation(old_window_start, id, removed);
        self.sorted.remove(slot);
        self.stats.phys_del += 1;
        if requeue {
            self.stats.risky += 1;
            self.mark_pending(id);
        }
    }

    /// One left rotation: window_start -> window_start + 1. Updates the factors and the sorted list.
    /// `best` = fewest runs found so far (used to decide whether the re-insertions can be skipped).
    fn rotate_left_once(&mut self, prune: bool, best: usize) {
        let (text_len, old_start) = (self.text_len, self.window_start);
        self.stats.steps += 1;
        let max_lcp = self.sorted.max_lcp; // read before any deletion changes it
        self.stats.lmax_sum += max_lcp as u64;
        let mut removed = std::mem::take(&mut self.removed);
        let mut pieces = std::mem::take(&mut self.pieces);
        let mut dead_owners = std::mem::take(&mut self.dead_owners);
        let mut born_owners = std::mem::take(&mut self.born_owners);
        removed.clear();
        pieces.clear();
        dead_owners.clear();
        born_owners.clear();

        // ---- A. new factorization: split the first factor, append the moved char, merge at the back ----
        let old_count = self.factors.len();
        let first = self.factors.pop_front().unwrap();
        // The first block survives iff it has a second copy: a merged factor starts with a letter
        // smaller than the moved char, the copy starts with that char (I1), so it is never merged.
        let first_block_survives = old_count > 1 && self.same_factor(first, self.factors[0]);
        // pieces = the Lyndon factors of first[1..], by walking the Lyndon table (I2)
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
        let mut back = Factor {
            start: old_start + text_len,
            len: 1,
        };
        self.factors.push_back(back);
        let mut old_left = old_count - 1; // old factors between the pieces and `back`
        while old_left > 0 {
            let before = self.factors[self.factors.len() - 2];
            if !factor_is_smaller(self.text, text_len, before, back) {
                break;
            }
            self.factors.pop_back();
            removed.push(before); // the merged factors are collected from the back: descending
            old_left -= 1;
            back = Factor {
                start: before.start,
                len: before.len + back.len,
            };

            *self.factors.back_mut().unwrap() = back;
        }
        // (I1: if no old factor is left, the factor before `back` is a piece, which is never smaller)
        debug_assert!(
            old_left > 0
                || self.factors.len() < 2
                || !factor_is_smaller(
                    self.text,
                    text_len,
                    self.factors[self.factors.len() - 2],
                    back
                )
        );
        removed.push(first);
        removed.reverse(); // ascending: first, then the merged factors in text order
        self.window_start = old_start + 1;

        // ---- B. which owners die, which are born (both lists ascending, disjoint) ----
        // Dead owners: the first factor without a second copy, and every merged block (by its last
        // copy). The first position of the window leaves at the front and re-enters at the back,
        // so the first factor's positions are [old_start+1, first.end) and old_start+text_len.
        let dead_first = !first_block_survives;
        if dead_first && first.len > 1 {
            dead_owners.push(Factor {
                start: old_start + 1,
                len: first.len - 1,
            });
        }
        // If `back` equals the old factor before it, no new factor is born: that block just grows,
        // and its entries (already in the list) are renamed to the new last copy, see below.
        let mut grown_block = None;
        if old_left > 0 {
            let before_back = self.factors[self.factors.len() - 2];
            if self.same_factor(before_back, back) {
                grown_block = Some(before_back);
                self.stats.growth += 1;
            }
        }
        for k in 1..removed.len() {
            if k + 1 == removed.len() || !self.same_factor(removed[k], removed[k + 1]) {
                dead_owners.push(removed[k]);
            }
        }
        if dead_first {
            dead_owners.push(Factor {
                start: old_start + text_len,
                len: 1,
            });
        }
        // Born owners: the last piece of each block of equal pieces, and `back` unless it only extends a block.
        for k in 0..pieces.len() {
            if k + 1 == pieces.len() || !self.same_factor(pieces[k], pieces[k + 1]) {
                born_owners.push(pieces[k]);
            }
        }
        if grown_block.is_none() {
            born_owners.push(back);
        }

        // ---- C. delete the entries that exist only before the rotation ----
        // (positions of dead owners that are not positions of born owners)
        let mut dead_count = 0u64;
        for d in &dead_owners {
            for_each_uncovered(d.start, d.end(), &born_owners, |q| {
                self.delete_entry(q % text_len, false, old_start, &removed);
                dead_count += 1;
            });
        }
        self.stats.dead += dead_count;

        // ---- D. entries that exist before and after ("kept") but may change place ----
        // A kept entry only moves if its word changes within max_lcp characters, i.e. it lies within
        // max_lcp of the end of its old or new factor; its letter changes at the first position of a
        // factor. Such entries are deleted now (old order) and inserted later (new order).
        let owners_contain = |owners: &[Factor], q: usize| {
            let k = owners.partition_point(|f| f.start <= q);
            k > 0 && q < owners[k - 1].end()
        };
        for f in dead_owners.iter().chain(born_owners.iter()) {
            let near_end = f.end().saturating_sub(max_lcp).max(f.start)..f.end();
            for q in near_end.chain([f.start]) {
                if owners_contain(&dead_owners, q) && owners_contain(&born_owners, q) {
                    self.delete_entry(q % text_len, true, old_start, &removed);
                }
            }
        }

        // ---- E. entries that exist only after the rotation: queue them for insertion ----
        for b in &born_owners {
            for_each_uncovered(b.start, b.end(), &dead_owners, |q| {
                self.mark_pending(q % text_len);
                self.stats.ins += 1;
            });
        }

        // ---- E2. a grown block: its entries move to the new last copy, ids shift by |back| ----
        // The slots are located first, in the unchanged list; only then are the ids changed (a lookup
        // must never see an entry whose id already points into the new structure). The order of the
        // entries does not change: they are the same conjugates of the same factor.
        if let Some(old_last) = grown_block {
            let mut renames = std::mem::take(&mut self.renames);
            renames.clear();
            for k in 0..back.len {
                let (old_id, new_id) =
                    ((old_last.start + k) % text_len, (back.start + k) % text_len);
                let slot = if self.is_pending[old_id] {
                    None
                } else {
                    Some(self.locate_entry_before_rotation(old_start, old_id, &removed))
                };
                renames.push((old_id, new_id, slot));
            }
            for &(old_id, new_id, slot) in &renames {
                match slot {
                    Some(s) => self.sorted.chunks[s.0][s.1].id = new_id as u32,
                    None => {
                        self.unmark_pending(old_id);
                        self.mark_pending(new_id);
                    }
                }
            }
            self.renames = renames;
        }

        self.removed = removed;
        self.pieces = pieces;
        self.dead_owners = dead_owners;
        self.born_owners = born_owners;

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
    let text_len = x.len();
    let (mut out, mut i) = (Vec::new(), 0);
    while i < text_len {
        let (mut j, mut k) = (i + 1, i);
        while j < text_len && x[k] <= x[j] {
            if x[k] < x[j] {
                k = i;
            } else {
                k += 1;
            }
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
    if v.is_empty() {
        0
    } else {
        1 + v.windows(2).filter(|w| w[0] != w[1]).count()
    }
}

impl<'a> RotationSweep<'a> {
    fn letters(&self) -> Vec<u8> {
        self.sorted
            .chunks
            .iter()
            .flatten()
            .map(|e| e.letter)
            .collect()
    }
    /// Checks that the sorted list is really sorted, that every stored letter/LCP/counter is right.
    fn check_invariants(&mut self) -> Result<(), String> {
        let (text_len, window_start, text) = (self.text_len, self.window_start, self.text);
        let flat: Vec<Entry> = self.sorted.chunks.iter().flatten().copied().collect();
        if flat.len() != self.sorted.len {
            return Err("len".into());
        }
        if self
            .sorted
            .chunks
            .iter()
            .any(|c| c.is_empty() || c.len() > 2 * CHUNK)
        {
            return Err("chunk size".into());
        }
        let mut scratch_stats = ChainStats::default();
        let mut max_lcp = 0;
        for j in 0..flat.len() {
            let conj = conjugate_of_entry(text_len, flat[j].id, &self.factors, window_start);
            if self.predecessor_letter(&conj) != flat[j].letter {
                return Err(format!("letter at {}", j));
            }
            if j + 1 < flat.len() {
                let next_conj =
                    conjugate_of_entry(text_len, flat[j + 1].id, &self.factors, window_start);
                let (lcp, order) =
                    compare_conjugates(text, text_len, &conj, &next_conj, &mut scratch_stats);
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
        if flat
            .windows(2)
            .filter(|w| w[0].letter != w[1].letter)
            .count()
            != self.sorted.letter_changes
        {
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
    if text_len == 0 {
        return (0, 0, ChainStats::default());
    }
    if text_len == 1 {
        return (0, 1, ChainStats::default());
    }
    let mut sweep = RotationSweep::new(text, circular_sa);
    let (mut best, mut best_shift) = (sweep.runs(), 0);
    for _ in 0..text_len - 1 {
        let best_before = best;
        sweep.rotate_left_once(prune, best);
        // the run count is exact only when nothing is pending
        if sweep.pending_count == 0 && sweep.runs() < best {
            best = sweep.runs();
            best_shift = sweep.window_start;
        }
        if verify {
            let x: Vec<u8> = (0..text_len)
                .map(|k| text[(sweep.window_start + k) % text_len])
                .collect();
            let expected = oracle_letters(&x);
            let true_runs = runs_of_letters(&expected);
            let mut flushed = sweep.clone();
            flushed.flush_pending();
            flushed.check_invariants().unwrap_or_else(|e| {
                panic!("invariant after step to {}: {}", sweep.window_start, e)
            });
            assert_eq!(
                flushed.letters(),
                expected,
                "letters differ at rotation {} (pending {})",
                sweep.window_start,
                sweep.pending_count
            );
            if sweep.pending_count > 0 {
                assert!(
                    sweep.runs() <= true_runs,
                    "partial list has more runs than the true count"
                );
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

/// Index (among the conjugates of the least rotation of the primitive root) of the first conjugate
/// with the fewest BBWT runs. `prune` enables the certified skipping.
pub fn optimal_bbwt_chain(input: &[u8], prune: bool) -> usize {
    let (start, period) = shiloach(input);
    let text: Vec<u8> = (0..period)
        .map(|k| input[(start + k) % input.len()])
        .collect();
    let circular_sa = circular_suffix_array(&text); // initial order of the text_len rotations: the only suffix-array computation of the run
    chain_min_runs_sa(&text, Some(&circular_sa), prune, false).0
}

/// Same sweep, cross-checked against the oracle at every step when n <= 400, with statistics.
pub fn online_bbwt_check(input: &[u8]) -> usize {
    let (start, period) = shiloach(input);
    let text: Vec<u8> = (0..period)
        .map(|k| input[(start + k) % input.len()])
        .collect();
    let verify = text.len() <= 400;
    let circular_sa = circular_suffix_array(&text);
    let (shift, runs, st) = chain_min_runs_sa(&text, Some(&circular_sa), true, verify);
    let steps = st.steps.max(1) as f64;
    println!("n={} verified={} best shift={} runs={} | per step: risky {:.1} deleted {:.1} reinserted {:.1} cmps {:.0} | deferred {} of {} steps, block growths {}",
        text.len(), verify, shift, runs, st.risky as f64 / steps, st.phys_del as f64 / steps, st.phys_ins as f64 / steps,
        st.cmp_calls as f64 / steps, st.deferred_steps, st.steps, st.growth);
    shift
}
