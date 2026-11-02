use rand::RngExt;

// generates a random binary string
pub fn generate_random_binary_string(length: usize) -> String {
    let mut rng = rand::rng();
    let mut result = String::with_capacity(length);
    for _ in 0..length {
        let bit = rng.random_range(0..=1);
        result.push_str(&bit.to_string());
    }
    result
}