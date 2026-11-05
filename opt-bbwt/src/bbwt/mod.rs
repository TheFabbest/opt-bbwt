pub mod bbwt;
pub mod naive;

pub fn bbwt(text: &[u8]) -> Vec<u8> {
    self::bbwt::bbwt(text)
}

pub fn bbwt_into(text: &[u8], output: &mut [u8]) {
    self::bbwt::bbwt_into(text, output)
}

pub fn unbbwt(text: &[u8]) -> Vec<u8> {
    self::bbwt::unbbwt(text)
}

pub fn unbbwt_into(text: &[u8], output: &mut [u8]) {
    self::bbwt::unbbwt_into(text, output)
}

pub mod variants {
    pub use super::bbwt::{bbwt, bbwt_into, unbbwt, unbbwt_into};
    pub use super::naive::naive_bbwt;
}