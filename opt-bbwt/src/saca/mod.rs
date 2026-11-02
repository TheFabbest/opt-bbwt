pub mod optimalsaca;
pub mod circularsais;
pub mod sais;

pub use self::optimalsaca::{optimal_suffix_array_bytes};
pub use self::circularsais::circular_suffix_array;
pub use self::sais::sais;