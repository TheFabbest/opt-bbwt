pub fn omega_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    // Compare their infinite expansions a.iter().cycle() and b.iter().cycle()
    let len_lcm = a.len() * b.len(); // or compare up to a.len() * b.len() steps
    let inf_a = a.iter().cycle().take(len_lcm);
    let inf_b = b.iter().cycle().take(len_lcm);
    inf_a.cmp(inf_b)
}
