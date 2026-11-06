use rand::RngExt;

/// Generates a random binary string of the requested length.
pub fn generate_random_binary_string(length: usize) -> String {
    let mut rng = rand::rng();
    let mut result = String::with_capacity(length);
    for _ in 0..length {
        let bit = rng.random_range(0..=1);
        result.push_str(&bit.to_string());
    }
    result
}
