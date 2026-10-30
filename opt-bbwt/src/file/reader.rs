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
    let path = format!("./corpora/{}", corpus_name);
    let bytes = read_file_to_bytes(&path)?;
    Ok(bytes[..length].to_vec())
}