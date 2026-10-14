use std::cmp::Ordering;
use std::fmt;

/// Little-endian 384-bit unsigned integer.  The Haskell implementation uses
/// Word320; the extra limb safely holds sums of four phase labels.
#[derive(Clone, Copy, Default, Eq, PartialEq, Hash)]
pub struct U384(pub [u64; 6]);

impl U384 {
    pub const ZERO: Self = Self([0; 6]);
    pub const ONE: Self = Self([1, 0, 0, 0, 0, 0]);

    pub fn from_u64(x: u64) -> Self {
        Self([x, 0, 0, 0, 0, 0])
    }

    pub fn pow2(bit: usize) -> Self {
        assert!(bit < 384);
        let mut out = Self::ZERO;
        out.0[bit / 64] = 1u64 << (bit % 64);
        out
    }

    pub fn is_zero(self) -> bool {
        self.0.iter().all(|&x| x == 0)
    }

    pub fn bit_len(self) -> usize {
        for i in (0..6).rev() {
            if self.0[i] != 0 {
                return i * 64 + (64 - self.0[i].leading_zeros() as usize);
            }
        }
        0
    }

    pub fn bit(self, i: usize) -> bool {
        ((self.0[i / 64] >> (i % 64)) & 1) != 0
    }

    pub fn set_bit(&mut self, i: usize) {
        self.0[i / 64] |= 1u64 << (i % 64);
    }

    pub fn checked_add(self, rhs: Self) -> Option<Self> {
        let mut out = [0u64; 6];
        let mut carry = 0u128;
        for (i, slot) in out.iter_mut().enumerate() {
            let z = self.0[i] as u128 + rhs.0[i] as u128 + carry;
            *slot = z as u64;
            carry = z >> 64;
        }
        (carry == 0).then_some(Self(out))
    }

    pub fn add(self, rhs: Self) -> Self {
        self.checked_add(rhs).expect("U384 addition overflow")
    }

    pub fn sub(self, rhs: Self) -> Self {
        assert!(self >= rhs, "U384 subtraction underflow");
        let mut out = [0u64; 6];
        let mut borrow = 0u128;
        for (i, slot) in out.iter_mut().enumerate() {
            let a = self.0[i] as u128;
            let b = rhs.0[i] as u128 + borrow;
            if a >= b {
                *slot = (a - b) as u64;
                borrow = 0;
            } else {
                *slot = ((1u128 << 64) + a - b) as u64;
                borrow = 1;
            }
        }
        debug_assert_eq!(borrow, 0);
        Self(out)
    }

    pub fn checked_mul_small(self, rhs: u64) -> Option<Self> {
        let mut out = [0u64; 6];
        let mut carry = 0u128;
        for (i, slot) in out.iter_mut().enumerate() {
            let z = self.0[i] as u128 * rhs as u128 + carry;
            *slot = z as u64;
            carry = z >> 64;
        }
        (carry == 0).then_some(Self(out))
    }

    pub fn mul_small(self, rhs: u64) -> Self {
        self.checked_mul_small(rhs).expect("U384 multiplication overflow")
    }

    pub fn checked_mul(self, rhs: Self) -> Option<Self> {
        let mut out = [0u64; 6];
        for i in 0..6 {
            let mut carry = 0u128;
            for j in 0..(6 - i) {
                let k = i + j;
                let z = out[k] as u128 + self.0[i] as u128 * rhs.0[j] as u128 + carry;
                out[k] = z as u64;
                carry = z >> 64;
            }
            if carry != 0 {
                return None;
            }
            for j in (6 - i)..6 {
                if self.0[i] != 0 && rhs.0[j] != 0 {
                    return None;
                }
            }
        }
        Some(Self(out))
    }

    pub fn mul(self, rhs: Self) -> Self {
        self.checked_mul(rhs).expect("U384 multiplication overflow")
    }

    pub fn div_rem_small(self, rhs: u64) -> (Self, u64) {
        assert_ne!(rhs, 0);
        let mut out = [0u64; 6];
        let mut rem = 0u128;
        for i in (0..6).rev() {
            let z = (rem << 64) | self.0[i] as u128;
            out[i] = (z / rhs as u128) as u64;
            rem = z % rhs as u128;
        }
        (Self(out), rem as u64)
    }

    pub fn div_small(self, rhs: u64) -> Self {
        self.div_rem_small(rhs).0
    }

    pub fn shl1(self) -> Self {
        let mut out = [0u64; 6];
        let mut carry = 0u64;
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = (self.0[i] << 1) | carry;
            carry = self.0[i] >> 63;
        }
        assert_eq!(carry, 0, "U384 shift overflow");
        Self(out)
    }

    pub fn div_rem(self, rhs: Self) -> (Self, Self) {
        assert!(!rhs.is_zero());
        if self < rhs {
            return (Self::ZERO, self);
        }
        let mut q = Self::ZERO;
        let mut r = Self::ZERO;
        for i in (0..self.bit_len()).rev() {
            r = r.shl1();
            if self.bit(i) {
                r.0[0] |= 1;
            }
            if r >= rhs {
                r = r.sub(rhs);
                q.set_bit(i);
            }
        }
        (q, r)
    }

    pub fn div(self, rhs: Self) -> Self {
        self.div_rem(rhs).0
    }

    pub fn modulo(self, modulus: Self) -> Self {
        self.div_rem(modulus).1
    }

    pub fn log2_f64(self) -> f64 {
        assert!(!self.is_zero());
        let bits = self.bit_len();
        if bits <= 53 {
            return (self.0[0] as f64).log2();
        }
        let shift = bits - 53;
        let top = self.shr_to_u64(shift);
        (top as f64).log2() + shift as f64
    }

    pub fn to_f64(self) -> f64 {
        if self.is_zero() {
            0.0
        } else {
            2f64.powf(self.log2_f64())
        }
    }

    pub fn to_u64(self) -> Option<u64> {
        self.0[1..].iter().all(|&x| x == 0).then_some(self.0[0])
    }

    fn shr_to_u64(self, shift: usize) -> u64 {
        let limb = shift / 64;
        let off = shift % 64;
        let mut v = self.0[limb] >> off;
        if off != 0 && limb + 1 < 6 {
            v |= self.0[limb + 1] << (64 - off);
        }
        v
    }

    pub fn from_decimal(s: &str) -> Result<Self, String> {
        let mut out = Self::ZERO;
        if s.is_empty() {
            return Err("empty integer".into());
        }
        for c in s.bytes() {
            if !c.is_ascii_digit() {
                return Err(format!("invalid decimal digit {:?}", c as char));
            }
            out = out
                .checked_mul_small(10)
                .and_then(|x| x.checked_add(Self::from_u64((c - b'0') as u64)))
                .ok_or_else(|| "decimal integer exceeds 384 bits".to_string())?;
        }
        Ok(out)
    }
}

impl Ord for U384 {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in (0..6).rev() {
            match self.0[i].cmp(&other.0[i]) {
                Ordering::Equal => {}
                x => return x,
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for U384 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for U384 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return f.write_str("0");
        }
        let mut n = *self;
        let mut digits = Vec::new();
        while !n.is_zero() {
            let (q, r) = n.div_rem_small(10);
            digits.push((b'0' + r as u8) as char);
            n = q;
        }
        for c in digits.iter().rev() {
            write!(f, "{c}")?;
        }
        Ok(())
    }
}

impl fmt::Debug for U384 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self}")
    }
}

#[cfg(test)]
mod tests {
    use super::U384;

    #[test]
    fn arithmetic_and_decimal() {
        let a = U384::from_decimal("123456789012345678901234567890").unwrap();
        let b = U384::from_decimal("98765432109876543210").unwrap();
        let p = a.mul(b);
        assert_eq!(p.div_rem(b), (a, U384::ZERO));
        assert_eq!(a.to_string(), "123456789012345678901234567890");
    }

    #[test]
    fn high_bit_division() {
        let a = U384::pow2(320).add(U384::from_u64(17));
        let b = U384::pow2(80).add(U384::from_u64(3));
        let (q, r) = a.div_rem(b);
        assert_eq!(q.mul(b).add(r), a);
        assert!(r < b);
    }
}
