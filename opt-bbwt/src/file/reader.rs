// reads a file into [u8] and returns it
pub fn read_file_to_bytes(path: &str) -> std::io::Result<Vec<u8>> {
    use std::fs::File;
    use std::io::Read;

    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    Ok(buffer)
}

// coprora loader
pub fn load_corpus(corpus_name: &str, length: usize) -> std::io::Result<Vec<u8>> {
    load_file_prefix(&format!("./corpora/{}", corpus_name), length)
}

/// Reads at most `length` bytes from a file without loading the complete file.
pub fn load_file_prefix(path: &str, length: usize) -> std::io::Result<Vec<u8>> {
    use std::fs::File;
    use std::io::{Read, Seek, SeekFrom};

    let mut file = File::open(path)?;
    let available = file.seek(SeekFrom::End(0))?;
    file.seek(SeekFrom::Start(0))?;

    let bytes_to_read = (length as u64).min(available);
    let capacity = usize::try_from(bytes_to_read).unwrap_or(usize::MAX);
    let mut buffer = Vec::with_capacity(capacity);
    file.take(bytes_to_read).read_to_end(&mut buffer)?;
    Ok(buffer)
}

#[cfg(test)]
mod tests {
    use super::load_file_prefix;

    #[test]
    fn prefix_is_bounded_by_requested_length() {
        let path = format!("{}/../corpora/dna", env!("CARGO_MANIFEST_DIR"));
        let bytes = load_file_prefix(&path, 8).expect("corpus should be available");
        assert_eq!(bytes.len(), 8);
    }

    #[test]
    fn zero_length_prefix_is_empty() {
        let path = format!("{}/../corpora/dna", env!("CARGO_MANIFEST_DIR"));
        assert!(load_file_prefix(&path, 0)
            .expect("corpus should be available")
            .is_empty());
    }
}
