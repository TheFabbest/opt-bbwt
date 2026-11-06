use std::io;
use std::path::{Path, PathBuf};

const PIZZA_CHILI_URL: &str = "https://pizzachili.dcc.uchile.cl/texts";
const CORPORA: &[(&str, &str)] = &[
    ("nlang", "english"),
    ("code", "sources"),
    ("dna", "dna"),
    ("music", "pitches"),
    ("protein", "proteins"),
    ("xml", "dblp.xml"),
];

/// Reads a complete file into a byte vector.
pub fn read_file_to_bytes(path: &str) -> std::io::Result<Vec<u8>> {
    use std::fs::File;
    use std::io::Read;

    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    Ok(buffer)
}

/// Loads at most `length` bytes from a cached Pizza&Chili corpus.
///
/// Missing corpora are downloaded and decompressed before being read.
pub fn load_corpus(corpus_name: &str, length: usize) -> std::io::Result<Vec<u8>> {
    let path = corpus_path(corpus_name)?;
    if !path.exists() {
        download_corpus(corpus_name)?;
    }
    load_file_prefix(path.to_str().unwrap(), length)
}

/// Downloads and caches a Pizza&Chili corpus under the shared corpus directory.
///
/// `corpus_name` must be one of `english`, `sources`, `dna`, `pitches`,
/// `proteins`, or `dblp.xml`.
pub fn download_corpus(corpus_name: &str) -> std::io::Result<PathBuf> {
    let (collection, canonical_name) = corpus_info(corpus_name)?;
    let path = corpus_root().join(canonical_name);
    if path.exists() {
        return Ok(path);
    }

    std::fs::create_dir_all(path.parent().unwrap())?;
    let url = format!("{PIZZA_CHILI_URL}/{collection}/{canonical_name}.gz");
    let response = ureq::get(&url)
        .call()
        .map_err(|error| io::Error::other(format!("failed to download {url}: {error}")))?;
    let temporary_path = path.with_extension("download");
    let result = (|| {
        let mut decoder = flate2::read::GzDecoder::new(response.into_reader());
        let mut output = std::fs::File::create(&temporary_path)?;
        io::copy(&mut decoder, &mut output)?;
        std::fs::rename(&temporary_path, &path)?;
        Ok(path.clone())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary_path);
    }
    result
}

fn corpus_path(corpus_name: &str) -> io::Result<PathBuf> {
    let (_, canonical_name) = corpus_info(corpus_name)?;
    Ok(corpus_root().join(canonical_name))
}

fn corpus_info(corpus_name: &str) -> io::Result<(&'static str, &'static str)> {
    CORPORA
        .iter()
        .find(|(_, name)| *name == corpus_name)
        .copied()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "unknown corpus name; see the README for supported names",
            )
        })
}

fn corpus_root() -> PathBuf {
    std::env::var_os("OPT_BBWT_CORPORA_DIR").map_or_else(
        || {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .join("corpora")
        },
        PathBuf::from,
    )
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
    use super::{corpus_path, load_file_prefix};
    use std::fs;

    #[test]
    fn prefix_is_bounded_by_requested_length() {
        let path = std::env::temp_dir().join("opt-bbwt-reader-test");
        fs::write(&path, b"abcdefgh").expect("test file should be writable");
        let bytes =
            load_file_prefix(path.to_str().unwrap(), 8).expect("test file should be available");
        fs::remove_file(path).expect("test file should be removable");
        assert_eq!(bytes.len(), 8);
    }

    #[test]
    fn zero_length_prefix_is_empty() {
        let path = std::env::temp_dir().join("opt-bbwt-reader-empty-test");
        fs::write(&path, b"abcdefgh").expect("test file should be writable");
        assert!(load_file_prefix(path.to_str().unwrap(), 0)
            .expect("test file should be available")
            .is_empty());
        fs::remove_file(path).expect("test file should be removable");
    }

    #[test]
    fn invalid_corpus_names_are_rejected() {
        assert!(corpus_path("../dna").is_err());
        assert!(corpus_path("dna.50MB").is_err());
    }

    #[test]
    fn corpus_names_resolve_to_full_files() {
        assert!(corpus_path("dna").unwrap().ends_with("dna"));
        assert!(corpus_path("sources").unwrap().ends_with("sources"));
        assert!(corpus_path("dblp.xml").unwrap().ends_with("dblp.xml"));
    }
}
