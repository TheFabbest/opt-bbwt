//! File and test-data helpers.

pub mod generator;
pub mod reader;

pub use generator::generate_random_binary_string;
pub use reader::{load_corpus, load_file_prefix, read_file_to_bytes};
