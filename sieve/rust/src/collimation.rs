use crate::uint::U384;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, VecDeque};

pub type Quad = [u32; 4];

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct Pair {
    value: U384,
    i: u32,
    j: u32,
}

struct AscPairs<'a> {
    a: &'a [U384],
    b: &'a [U384],
    heap: BinaryHeap<Reverse<Pair>>,
}

impl<'a> AscPairs<'a> {
    fn new(a: &'a [U384], b: &'a [U384]) -> Self {
        let mut heap = BinaryHeap::new();
        if !b.is_empty() {
            for (i, &x) in a.iter().enumerate() {
                heap.push(Reverse(Pair { value: x.add(b[0]), i: i as u32, j: 0 }));
            }
        }
        Self { a, b, heap }
    }
}

impl Iterator for AscPairs<'_> {
    type Item = Pair;
    fn next(&mut self) -> Option<Pair> {
        let Reverse(p) = self.heap.pop()?;
        let next_j = p.j as usize + 1;
        if next_j < self.b.len() {
            self.heap.push(Reverse(Pair {
                value: self.a[p.i as usize].add(self.b[next_j]),
                i: p.i,
                j: next_j as u32,
            }));
        }
        Some(p)
    }
}

struct DescPairs<'a> {
    a: &'a [U384],
    b: &'a [U384],
    heap: BinaryHeap<Pair>,
}

impl<'a> DescPairs<'a> {
    fn new(a: &'a [U384], b: &'a [U384]) -> Self {
        let mut heap = BinaryHeap::new();
        if let Some(&last) = b.last() {
            for (i, &x) in a.iter().enumerate() {
                heap.push(Pair { value: x.add(last), i: i as u32, j: (b.len() - 1) as u32 });
            }
        }
        Self { a, b, heap }
    }
}

impl Iterator for DescPairs<'_> {
    type Item = Pair;
    fn next(&mut self) -> Option<Pair> {
        let p = self.heap.pop()?;
        if p.j > 0 {
            let j = p.j - 1;
            self.heap.push(Pair {
                value: self.a[p.i as usize].add(self.b[j as usize]),
                i: p.i,
                j,
            });
        }
        Some(p)
    }
}

fn emit_four<F: FnMut(U384, Quad)>(
    interval: U384,
    quotient: U384,
    vectors: [&[U384]; 4],
    mut emit: F,
) {
    let start = quotient.mul(interval);
    let end = start.add(interval);
    let mut left = AscPairs::new(vectors[0], vectors[1]).peekable();
    let mut right = DescPairs::new(vectors[2], vectors[3]).peekable();
    let mut window: VecDeque<Pair> = VecDeque::new();

    while let Some(first) = left.next() {
        if first.value >= end {
            break;
        }
        let value = first.value;
        let mut group = vec![first];
        while left.peek().is_some_and(|p| p.value == value) {
            group.push(left.next().unwrap());
        }

        let lower = if start > value { start.sub(value) } else { U384::ZERO };
        let upper = end.sub(U384::ONE).sub(value);
        while window.front().is_some_and(|p| p.value > upper) {
            window.pop_front();
        }
        while right.peek().is_some_and(|p| p.value > upper) {
            right.next();
        }
        while right.peek().is_some_and(|p| p.value >= lower) {
            window.push_back(right.next().unwrap());
        }
        for x in &group {
            for y in &window {
                let total = x.value.add(y.value);
                emit(total.sub(start), [x.i, x.j, y.i, y.j]);
            }
        }
    }
}

pub fn collimate2(interval: U384, quotient: U384, a: &[U384], b: &[U384]) -> Vec<U384> {
    let start = quotient.mul(interval);
    let end = start.add(interval);
    let mut out = Vec::new();
    for &x in a {
        if x >= end {
            break;
        }
        let lower = if start > x {
            start.sub(x)
        } else {
            U384::ZERO
        };
        let upper = end.sub(x);
        let lo = b.partition_point(|&y| y < lower);
        let hi = b.partition_point(|&y| y < upper);
        for &y in &b[lo..hi] {
            out.push(x.add(y).sub(start));
        }
    }
    out.sort_unstable();
    out
}

pub fn collimate4(interval: U384, quotient: U384, vectors: [&[U384]; 4]) -> Vec<U384> {
    let mut out = Vec::new();
    emit_four(interval, quotient, vectors, |phase, _| out.push(phase));
    out.sort_unstable();
    out
}

pub struct FourListMap {
    pub phases: Vec<U384>,
    pub inv: Vec<Quad>,
    pub pi: HashMap<Quad, usize>,
}

pub fn collimate4_mapped(interval: U384, quotient: U384, vectors: [&[U384]; 4]) -> FourListMap {
    let mut entries = Vec::new();
    emit_four(interval, quotient, vectors, |phase, quad| entries.push((phase, quad)));
    entries.sort_unstable();
    let mut phases = Vec::with_capacity(entries.len());
    let mut inv = Vec::with_capacity(entries.len());
    let mut pi = HashMap::with_capacity(entries.len());
    for (j, (phase, quad)) in entries.into_iter().enumerate() {
        phases.push(phase);
        inv.push(quad);
        pi.insert(quad, j);
    }
    FourListMap { phases, inv, pi }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn u(xs: &[u64]) -> Vec<U384> { xs.iter().copied().map(U384::from_u64).collect() }

    #[test]
    fn binary_underflow_regression() {
        assert_eq!(collimate2(U384::from_u64(8), U384::ONE, &u(&[9]), &u(&[0])), u(&[1]));
    }

    #[test]
    fn binary_matches_exhaustive_with_duplicates() {
        let a = u(&[0, 1, 1, 4, 8]);
        let b = u(&[0, 2, 2, 5, 9]);
        for q in 0..5 {
            let mut expected = Vec::new();
            for &x in &a {
                for &y in &b {
                    let sum = x.add(y);
                    let start = U384::from_u64(q * 4);
                    let end = start.add(U384::from_u64(4));
                    if sum >= start && sum < end {
                        expected.push(sum.sub(start));
                    }
                }
            }
            expected.sort_unstable();
            assert_eq!(collimate2(U384::from_u64(4), U384::from_u64(q), &a, &b), expected);
        }
    }

    #[test]
    fn pair_streams_preserve_duplicates() {
        let a = u(&[0, 2, 2]);
        let b = u(&[1, 3]);
        let asc: Vec<_> = AscPairs::new(&a, &b).map(|x| x.value).collect();
        let desc: Vec<_> = DescPairs::new(&a, &b).map(|x| x.value).collect();
        assert_eq!(asc, u(&[1, 3, 3, 3, 5, 5]));
        assert_eq!(desc, u(&[5, 5, 3, 3, 3, 1]));
    }

    fn brute(interval: u64, q: u64, vs: [&[U384]; 4]) -> Vec<(U384, Quad)> {
        let start = q * interval;
        let end = start + interval;
        let mut out = Vec::new();
        for (i, &a) in vs[0].iter().enumerate() {
            for (j, &b) in vs[1].iter().enumerate() {
                for (k, &c) in vs[2].iter().enumerate() {
                    for (l, &d) in vs[3].iter().enumerate() {
                        let sum = a.add(b).add(c).add(d);
                        if sum >= U384::from_u64(start) && sum < U384::from_u64(end) {
                            out.push((sum.sub(U384::from_u64(start)), [i as u32,j as u32,k as u32,l as u32]));
                        }
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }

    #[test]
    fn mapped_matches_exhaustive() {
        let cases = [
            (4, 0, u(&[0,1,1]), u(&[0,2]), u(&[0,1]), u(&[0,2])),
            (4, 1, u(&[0,1,1]), u(&[0,2]), u(&[0,1]), u(&[0,2])),
            (3, 1, u(&[0,1,2]), u(&[0,1]), u(&[0,1,2]), u(&[0,1])),
            (5, 2, u(&[1,4,4]), u(&[0,3]), u(&[2,2,5]), u(&[1,4])),
        ];
        for (s, q, a, b, c, d) in cases {
            let vs = [&a[..], &b[..], &c[..], &d[..]];
            let expected = brute(s, q, vs);
            let got = collimate4_mapped(U384::from_u64(s), U384::from_u64(q), vs);
            let pairs: Vec<_> = got.phases.iter().copied().zip(got.inv.iter().copied()).collect();
            assert_eq!(pairs, expected);
            for (j, quad) in got.inv.iter().enumerate() {
                assert_eq!(got.pi[quad], j);
            }
        }
    }
}
