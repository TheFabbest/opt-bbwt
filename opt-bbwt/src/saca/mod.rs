pub mod optimalsaca;
pub mod circularsais;

pub use self::optimalsaca::{optimal_suffix_array_bytes, optimal_suffix_array_recursive};
pub use self::circularsais::circular_suffix_array;