pub mod duval;
pub mod shiloach;
pub mod conjugatesfactors;
pub mod lyndontables;
pub mod lcp;

pub use self::duval::duval;
pub use self::shiloach::Shiloach;
pub use self::conjugatesfactors::print_all_factors;
pub use self::lyndontables::{lyndon_table, lyndon_suffix_table, compute_lyn_optimized};
pub use self::lcp::lcp;

pub trait CanonicalRotationFinder {
    /// Finds the canonical (minimal) rotation of a slice.
    /// Returns the starting index and the length of the canonical rotation.
    /// index + length must be considered modulo the length of the slice.
    fn find<T: Ord>(&self, s: &[T]) -> (usize, usize);
}