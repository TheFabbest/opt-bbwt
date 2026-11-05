# onlinebbwt.rs — what every piece does

Line numbers refer to `onlinebbwt.rs`. `text_len` is the length of the text; `text` is the base string, a Lyndon word.

## The idea in five sentences
1. The BBWT of a string x has as many runs as the sorted list of "entries" has letter changes. An entry is one conjugate of one *distinct* Lyndon factor of x; its letter is the character just before it inside the factor.
2. We visit every conjugate of `text` by rotating left one character at a time: x = a·y becomes y·a.
3. One rotation changes only a few Lyndon factors: the first factor splits into pieces, and the moved character merges with factors at the back.
4. So instead of rebuilding the sorted list, we delete the entries of the factors that died, insert the entries of the new factors, and fix only the surviving entries whose position in the order could change.
5. If, after the deletions, the list already has at least `best` runs, this conjugate cannot win, so we postpone the insertions.

## Where the circular suffix array is computed
Only once per run: line 707 (`optimal_bbwt_chain`) and line 716 (`online_bbwt_check`), `circular_suffix_array(&text)`. `RotationSweep::new` (line 354) uses it as the initial sorted order of the entries. If you pass `None` it sorts the rotations by comparison instead (line 368). After that the order is only ever edited.

## Basic building blocks

**`mul_mod` (51)** — returns (a·b) mod (2^61 − 1). It multiplies in 128 bits, then adds the high bits to the low bits. Used only for hashing.

**`Factor` (60)** — one occurrence of a Lyndon factor: `text[start .. start+len]`. `start` is an *absolute* position: it can be larger than `text_len`, because the current string wraps around the end of `text`.

**`wrap_index(x, text_len)` (67)** — returns x mod text_len, assuming x < 3·text_len. Absolute positions can run past the end of `text`; this maps them back to a real index of `text`. It uses two comparisons instead of `%` because it runs in the innermost loop.

**`Conjugate` (74)** — "the factor at `factor_start`, `factor_len` long, rotated left by `offset`". This is the string an entry stands for.

**`compare_conjugates` (88)** — compares two conjugates and returns (length of their common prefix, Less/Equal/Greater).
How it works: each conjugate is treated as an infinite word that repeats its own factor forever. The loop reads one character from each side (line 96) and stops at the first difference. When one side reaches the end of its factor, it jumps back to the factor's start (lines 98–99). If the two words agree for `factor_len_a + factor_len_b` characters they are equal, so the loop never needs to go further (line 90).

**`conjugate_of_entry` (107)** — turns an entry id (a text position) into its `Conjugate`: computes the absolute position of the id inside the current window, then calls `conjugate_at`.

**`conjugate_at` (114)** — finds which factor contains an absolute position (binary search over the sorted factor list) and returns the conjugate that starts there.

**`conjugate_of_entry_before_rotation` (123)** — same as `conjugate_of_entry`, but for the window as it was *before* the current rotation. The factor list has already been updated, so it rebuilds the old factors from `removed`: if the position lies inside a removed factor, that factor is the old one; otherwise the position was untouched and the new list gives the right answer.

**`factor_is_smaller(f, g)` (138)** — returns true if f < g as plain strings in lexicographic order (a proper prefix counts as smaller). It compares character by character up to the shorter length; the first difference decides. It is used for the Lyndon merge: if f < g then f·g is a Lyndon word, so two neighbouring factors f, g with f < g are not a valid factorization and must be merged.

**`range_difference`, `range_intersection`, `range_contains` (149, 164, 174)** — set operations on sets of positions stored as sorted lists of half-open ranges (lo, hi). Difference = in a but not b. Intersection = in both. Contains = is position p in some range.

## The sorted entry list: `SortedEntries` (194)

It stores all entries in omega-sorted order, split into chunks of about 512 so that inserting or deleting moves at most ~1000 entries instead of all of them.

Each `Entry` (183) has: `id` (the text position that identifies it), `letter`, and `lcp_next` (longest common prefix with the *next* entry in the list).

It also keeps two counters up to date at all times:
- `letter_changes` = number of neighbouring pairs with different letters, so `runs() = letter_changes + 1` (206);
- `lcp_histogram` and `max_lcp` = how many neighbouring pairs have each LCP, and the largest one. `max_lcp` is the longest prefix any two entries can share.

| Function | What it does |
|---|---|
| `from_sorted` (219) | builds the list from entries that are already sorted, with their LCPs |
| `seek` (232) | binary search in two levels (over the first entry of each chunk, then inside one chunk). Returns the first slot whose entry is not smaller than the target. The caller supplies `compare_to_target(id)` |
| `neighbours_of_gap` (249) | the slots just before and just after a given slot, also across chunk borders |
| `remove` (259) | deletes an entry and repairs the counters. Key point: because the list is sorted, the LCP between the two new neighbours is simply `min(lcp(prev,e), lcp(e,next))`, so no string comparison is needed (line 270) |
| `insert` (292) | inserts an entry; the caller already knows its LCP with both neighbours. Splits the chunk if it grew over 1024 |

## The sweep state: `RotationSweep` (336)

| Field | Meaning |
|---|---|
| `text`, `text_len` | the base string and its length |
| `window_start` | the current conjugate is `text[window_start..] + text[..window_start]` |
| `prefix_hash`, `base_power` | hashes used to test whether two factors are the same string |
| `factors` | the Lyndon factors of the current conjugate, left to right |
| `sorted` | the sorted entries |
| `pending`, `is_pending`, `pending_count` | entries that must be inserted but were postponed |
| `stats` | counters for reporting only |

Methods:
- **`new` (354)** — builds the starting state: `text` is Lyndon, so it is one factor and its entries are its `text_len` rotations. Takes their order from the circular suffix array, computes neighbour LCPs with `compare_conjugates`, and builds `sorted`.
- **`factor_hash` / `same_factor` (381, 385)** — two factors are treated as equal if they have the same length and hash. There is no character check.
- **`predecessor_letter` (389)** — the BBWT letter of a conjugate: the character just before it inside its factor (for offset 0, the factor's last character).
- **`mark_pending` / `unmark_pending` (392, 395)** — add or remove an id from the postponed set.
- **`insert_entry` (401)** — inserts one entry into `sorted`: binary search with `seek` for its place, then measure its LCP with both neighbours, then `insert`.
- **`flush_pending` (419)** — inserts every pending id, skipping ids that were cancelled in the meantime.
- **`lyndon_factors_of(start, len)` (433)** — Duval's algorithm on `text[start..start+len]`: returns its Lyndon factors. (`i` = start of the group being built, `j` = scan position, `k` = position `j` is compared with.)
- **`locate_entry_before_rotation` (451)** — finds where an entry sits in the sorted list as it was before this rotation, by binary search. It answers "Equal" for the entry itself, so the search stops there and never compares an entry with itself (which would cost a whole factor length).

## `rotate_left_once` (464): the whole algorithm for one rotation

**A. New factorization (468).**
1. Pop the first factor (`first`). It loses its first character.
2. If it has more than one character, run Duval on the rest (`lyndon_factors_of`, line 474); the pieces are pushed at the front and recorded as `added`.
3. The moved character becomes a new one-character factor at the back (`moved_char`).
4. While the second-to-last factor is smaller than the last (`factor_is_smaller`), merge them (lines 483–493). A factor created in this same step and merged away again cancels out of `added`; an old one goes to `removed`.
5. `first_block_survives` (496): the first factor's value still exists if it had an equal second copy and that copy was not merged away.

**B. Which distinct factors die or are born (501).**
- `dead_factors`: for each block of equal removed factors, the last occurrence. Equal factors form a block; its entries live in the last occurrence. A block dies completely unless it is the surviving first block.
- `growth` (512): if the new back factor equals the factor before it, the block just grows. Nothing is born; the block's entries must move to the new last occurrence.
- `born_factors` (522): for each block of equal added factors, the last occurrence.

**C. Positions of entries that die, are born, or stay (534).** Each factor owns the positions `start .. start+len` (its entries). Position `old_start` leaves the window at the front and re-enters at the back as `old_start + text_len` (line 537). Then:
- `dead_ranges` = old ranges not in new ranges: delete for good;
- `born_ranges` = new ranges not in old ranges: insert;
- `kept_ranges` = in both: the entry stays but its conjugate changed, so its place in the order may change.

**D. Which kept entries can really move (548).** An entry's word changes only after its factor ends. If that happens after more than `max_lcp` characters, no other entry shares that long a prefix with it, so its place in the order cannot change. So only these are "risky" (lines 552–562): kept entries within `max_lcp` positions of the end of an old or new factor, and the first position of each such factor (its letter changes there).

**E. Delete and queue (567).**
- For each dead id and each risky id: if it is currently pending, it is not in the list, so skip (cancel the pending mark if it dies); otherwise find it with `locate_entry_before_rotation`, `remove` it, and, if risky, mark it pending again (lines 570–576).
- Born ids are marked pending (577).
- Growth (578): the ids of the old last occurrence are renamed to the new one. All locations are found first, then renamed, so no lookup sees a half-renamed list.

**F. Insert or skip (596).** After step E the list is missing some entries, so it is a subsequence of the true letter sequence, and a subsequence never has more runs. If `prune` is on and it already has `>= best` runs, this conjugate cannot beat the best, so the insertions stay pending. Otherwise `flush_pending` makes the list exact.

## Verification code (603–660)
`oracle_letters` rebuilds the BBWT letters of any string from scratch (slow, independent). `check_invariants` verifies that the list is sorted, that every stored letter and LCP is correct, and that the counters match. They run only when `verify` is true.

## Entry points (666–723)
- **`chain_min_runs_sa` (671)** — the driver. Builds the sweep (672), then calls `rotate_left_once` `text_len − 1` times (677). It records a new best only when `pending_count == 0`, because then the run count is exact (679). With `verify`, it compares with the oracle after every step.
- **`chain_min_runs` (666)** — same without a suffix array.
- **`optimal_bbwt_chain` (704)** — finds the least rotation and primitive root with `shiloach`, computes the circular suffix array (707), runs the sweep, returns the best shift.
- **`online_bbwt_check` (712)** — the same, verified against the oracle when `text_len <= 400`, and it prints statistics.