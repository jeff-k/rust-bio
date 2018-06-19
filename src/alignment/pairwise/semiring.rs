use std::ops::{Add, Mul};
use std::cmp::min;
use std::fmt::Debug;

extern crate geo;
use geo::{Polygon, Point};
use geo::algorithm::convexhull;

use Tropical::*;

trait Semiring: Mul<Output=Self> + Add<Output=Self> + Sized {
    // additive identity, multiplicative absorptive
    fn identity() -> Self;
}

// Tropical semiring
impl Semiring for Tropical {
    // a ⊗ 0 = a
    fn identity() -> Tropical {
        N(0)
    }
}

#[derive(Debug, Copy, Clone)]
enum Tropical {
    Inf,
    N(u32),
}

impl Add for Tropical {
    type Output = Tropical;
    fn add(self, other: Tropical) -> Tropical {
        match (self, other) {
            (a, Inf) => a,
            (Inf, b) => b,
            (N(a), N(b)) => N(min(a, b)),
        }
    }
}

impl Mul for Tropical {
    type Output = Tropical;
    fn mul(self, other: Tropical) -> Tropical {
        match (self, other) {
            (_, Inf) => Inf,
            (Inf, _) => Inf,
            (N(a), N(b)) => N(a + b),
        }
    }
}

// Viterbi semiring for hidden markov models
struct Viterbi {
    p: f64,
}

// Polytope semiring for parametric alignment
// nb. in one dimension this is the tropical semiring

fn union(a: Polygon<f64>, b: Polygon<f64>) {
    // union of polygons
    Polygon::new(vec![], vec![]);
}

fn minkowski(a: Polygon<f64>, b: Polygon<f64>) {
    // Minkowski sum (dilation)
}

impl Add for Polygon<f64> {
    type Output = Polygon<f64>;
    fn add(self, other: Polygon<f64>) -> Polygon<f64> {
        // convex hull of union
        convexhull(union(self, other))
    }
}

impl Mul for Polygon<f64> {
    type Output = Polygon<f64>;
    fn mul(self, other: Polygon<f64>) -> Polygon<f64> {
        // Minkowski sum (dilation)
        let ps = Vec::new();
        for (p1, p2) in zip(self, other) {
            ps.push(p1 + p2); 
        }
        Polygon::new(ps, vec![]);
    }
}

impl Semiring for Polygon<f64> {
    fn identity() -> Polygon<f64> { Polygon::new(vec![], vec![]) }
}

#[derive(Debug)]
enum Op<C: Eq> {
    M(C, C), // match or substitution
    I, // insertion
    D, // deletion
    Dx, Ix, // gap open
}

fn levenshtein(op: Op<u8>) -> Tropical {
    match op {
        Op::D => N(1),
        Op::I => N(1),
        Op::Dx => N(1), // gap extension costs the same as gap open
        Op::Ix => N(1),
        Op::M(a, b) => { if a != b { N(1) } else { N(0) } },
    }
}

fn polyprop() -> Polytope {
    // construct alignment polytope
    let (gap, mismatch) = (1, 1);

    let score = |op: Op<u8>| -> Tropical {
        match op {
            Op::D => N(gap),
            Op::I => N(gap),
            Op::Dx => N(gap),
            Op::Ix => N(gap),
            Op::M(a, b) => { if a != b { N(mismatch) } else { N(0) }},
        }
    };
    align(&[1,2,3,4], &[1,2,3,4], score);
    Polygon::new()
}

fn align<T: Semiring + Copy + Clone + Debug>
                    (seq1: &[i32], seq2: &[i32], score: fn(Op<i32>) -> T) {
    let mut row: Vec<Vec<T>> = vec![vec![T::identity()];seq2.len()];
    row[0][0] = score(Op::M(seq1[0], seq2[0]));

    // defining the boundary condition and recursive structure, which is
    // further generalized with generalized algebraic dynamic proramming
    for j in 1..seq2.len() {
        for i in 1..seq1.len() {
            row[i][j]
                = (score(Op::M(seq1[i], seq2[j])) * row[i-1][j-1])
                + (score(Op::I) * row[i-1][j])
                + (score(Op::D) * row[j][i-1]);
        }
    }
    row[seq1.len()][row2.len()];
}
