/// Copies a periodic slice into a new vector.
pub fn from_period_slice(s: &[u8], start: usize, period_length: usize) -> Vec<u8> {
    let mut output = vec![0; period_length];
    output.resize(period_length, 0);
    from_period_slice_into(s, start, period_length, &mut output);
    output
}

/// Copies a periodic slice into a caller-provided buffer.
///
/// `output` must have length `period_length`.
///
/// # Panics
///
/// Panics if `output.len() != period_length` or if a non-empty output is
/// requested from an empty input.
pub fn from_period_slice_into(s: &[u8], start: usize, period_length: usize, output: &mut [u8]) {
    assert_eq!(
        output.len(),
        period_length,
        "periodic slice output buffer must match requested length"
    );
    if output.is_empty() {
        return;
    }
    assert!(
        !s.is_empty(),
        "cannot take a periodic slice from an empty input"
    );
    let start = start % s.len();
    for (offset, byte) in output.iter_mut().enumerate() {
        *byte = s[(start + offset) % s.len()];
    }
}

/// Returns a rotation of `s` starting at `start`.
pub fn rotate_string(s: &[u8], start: usize) -> Vec<u8> {
    let mut output = vec![0; s.len()];
    rotate_string_into(s, start, &mut output);
    output
}

/// Rotates a string into a caller-provided buffer.
///
/// `output` must have the same length as `s`.
///
/// # Panics
///
/// Panics if `output.len() != s.len()`.
pub fn rotate_string_into(s: &[u8], start: usize, output: &mut [u8]) {
    assert_eq!(
        output.len(),
        s.len(),
        "rotation output buffer must match input length"
    );
    if s.is_empty() {
        return;
    }
    let start = start % s.len();
    for (offset, byte) in output.iter_mut().enumerate() {
        *byte = s[(start + offset) % s.len()];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_period_slice_empty() {
        let s = vec![];
        let result = from_period_slice(&s, 0, 0);
        assert_eq!(result, vec![]);
    }

    #[test]
    fn test_from_period_slice_single_char() {
        let s = vec![1];
        let result = from_period_slice(&s, 0, 1);
        assert_eq!(result, vec![1]);
    }

    #[test]
    fn test_from_period_slice_multiple_chars() {
        let s = vec![1, 2, 3, 4, 5];
        let result = from_period_slice(&s, 4, 5);
        assert_eq!(result, vec![5, 1, 2, 3, 4]);
    }

    #[test]
    fn test_from_period_slice_periodic() {
        let s = vec![1, 2, 3, 4, 1, 2, 3, 4];
        let result = from_period_slice(&s, 5, 4);
        assert_eq!(result, vec![2, 3, 4, 1]);
    }

    #[test]
    fn test_from_period_slice_into() {
        let mut output = [0; 4];
        from_period_slice_into(&[1, 2, 3, 4], 3, 4, &mut output);
        assert_eq!(output, [4, 1, 2, 3]);
    }

    #[test]
    fn test_rotate_string_into() {
        let mut output = [0; 4];
        rotate_string_into(&[1, 2, 3, 4], 2, &mut output);
        assert_eq!(output, [3, 4, 1, 2]);
    }
}
