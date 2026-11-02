pub mod reader;
pub mod generator;

pub use reader::{read_file_to_bytes, load_corpus};
pub use generator::generate_random_binary_string;