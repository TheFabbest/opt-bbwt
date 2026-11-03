pub mod reader;
pub mod generator;

pub use reader::{load_corpus, load_file_prefix, read_file_to_bytes};
pub use generator::generate_random_binary_string;