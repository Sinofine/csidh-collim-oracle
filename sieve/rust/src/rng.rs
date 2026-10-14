// SplitMix64, written by Sebastiano Vigna.
// Adapted from the public-domain reference implementation:
// https://prng.di.unimi.it/splitmix64.c
//
// Used only as a deterministic PRNG for reproducible simulation;
// it is not cryptographically secure.

use crate::uint::U384;

#[derive(Clone)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    pub fn index(&mut self, n: usize) -> usize {
        assert!(n > 0);
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }

    pub fn unit_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    pub fn below(&mut self, modulus: U384) -> U384 {
        assert!(!modulus.is_zero());
        let bits = modulus.bit_len();
        loop {
            let mut limbs = [0u64; 6];
            for x in &mut limbs {
                *x = self.next_u64();
            }
            let used = (bits + 63) / 64;
            for x in &mut limbs[used..] {
                *x = 0;
            }
            let rem = bits % 64;
            if rem != 0 {
                limbs[used - 1] &= (1u64 << rem) - 1;
            }
            let x = U384(limbs);
            if x < modulus {
                return x;
            }
        }
    }
}
