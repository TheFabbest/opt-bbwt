pub fn lyndon_table(s: &[u8], rank: &[usize]) -> Vec<usize> {
    let len = s.len();
    let mut lyndon_table = vec![0; len];
    
    for i in (0..len).rev() {
        lyndon_table[i] = 1;
        let mut j = i + 1;
        while j < len && rank[i] < rank[j] {
            lyndon_table[i] += lyndon_table[j];
            j += lyndon_table[j];
        }
    }

    lyndon_table
}

/// Computes the Lyndon Table (Lyn) in O(n) utilizing the newRank shortcut 
/// to minimize expensive cache misses on the rank array.
pub fn compute_lyn_optimized(s: &[u8], rank: &[usize], lyn_s: &[usize]) -> Vec<usize> {
    let n = s.len();
    // 2. Build the Lyn table backwards using the newRank shortcut
    let mut lyn = vec![1; n];
    for i in (0..n).rev() {
        lyn[i] = 1;
        let mut j = i + 1;
        
        while j < n {
            // The O(1) mathematical shortcut logic:
            // If LynS[j-1] < j-i, we know rank[i] < rank[j] is FALSE.
            let can_merge = if lyn_s[j - 1] < (j - i) {
                false 
            } else {
                // Fallback to the actual rank comparison only when necessary
                rank[i] < rank[j]
            };

            if can_merge {
                lyn[i] += lyn[j];
                j += lyn[j];
            } else {
                break;
            }
        }
    }
    
    lyn
}

pub fn lyndon_suffix_table(s: &[u8]) -> Vec<usize> {
    let mut lyndon_suffix_table = vec![0; s.len()];
    lyndon_suffix_table[0] = 1;
    let mut per = 1;
    let mut i = 0;
    for j in 1..s.len() {
        if s[j] != s[i] {
            lyndon_suffix_table[j] = j + 1;
            per = j + 1;
            i = 0;
        } else {
            lyndon_suffix_table[j] = lyndon_suffix_table[i];
            i = (i + 1) % per;
        }
    }
    lyndon_suffix_table
}
