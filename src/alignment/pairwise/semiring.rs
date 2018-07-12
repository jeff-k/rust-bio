use std::ops::{Add, Mul};
use std::cmp::{max, Ordering, min};
use std::fmt::Debug;

use utils::TextSlice;

extern crate geo;
use geo::{Polygon, Point};
use geo::algorithm::convexhull;

use Tropical::*;

trait Semiring: Mul<Output=Self> + Add<Output=Self> + Sized {
    // additive identity, multiplicative absorbing
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
// nb. in one dimension this is the tropical semiring (well, it should be)

fn union(a: Polygon<f64>, b: Polygon<f64>) {
    // union of polygons
    Polygon::new(vec![], vec![]);
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

// polytope propagation
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
                    seq1: TextSlice, seq2: TextSlice, score: fn(Op<i32>) -> T {
    let mut row: Vec<Vec<T>> = vec![vec![T::identity()];seq2.len()];
    row[0][0] = score(Op::M(seq1[0], seq2[0]));
    let (m, n) = (seq1.len(), seq2.len());

    // initialize matrix
    let mut traceback: Vec<Vec<Cell>> =
            vec![vec![T::new(); n + 1]; m + 1];
    let mut ops: Vec<Op> = vec![];

    for i in 1..(m + 1) {
        traceback[i][0] = Cell { score: -1 * i as i32, op: Op::Del(None) };
    }
    for j in 1..(n + 1) {
        traceback[0][j] = Cell { score: -1 * j as i32, op: Op::Ins(None) };
    }

    traceback[0][0] = Cell { score: 0, op: Op::Match(None) };

    // store the last visited node in topological order so that
    // we can index into the end of the alignment
    for (i_p, r) in seq1.iter().enumerate() {
        let i = i_p + 1;
        // query base and index (traceback matrix rows)
        for (j_p, q) in seq2.iter().enumerate() {
            let j = j_p + 1;

            // match and deletion scores for first reference base
            let (mat, del) = if prevs.len() == 0 {
                (Cell { score: traceback[0][j - 1].score + (self.scoring)(*r, *q),
                        op: Op::Match(None) },
                 Cell { score: traceback[0][j].score - 1i32, op: Op::D })
            } else {
                let mut mat_max = Cell { score: T::Inf, op: Op::M };
                let mut del_max = Cell { score: T::Inf, op: Op::D };
                for prev_n in 0..prevs.len() {
                    let i_p: usize = prevs[prev_n].index() + 1;
                    mat_max = max(mat_max,
                        Cell { score: traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                op: Op::Match(Some((i_p - 1, i - 1)))});
                    del_max = max(del_max,
                        Cell { score: traceback[i_p][j].score - 1i32,
                               op: Op::Del(Some((i_p - 1, i)))});
                }
                (mat_max, del_max)
            };
            let score = max(mat, max(del, Cell { score: traceback[i][j - 1].score - 1i32, 
                                                 op: Op::Ins(Some(i - 1)) }));
            traceback[i][j] = score;
        }
    }
    
    // Now backtrack through the matrix to construct the optimal path
    while i > 0 && j > 0 {
        // push operation and edge corresponding to (one of the) optimal
        // routes
        //println!("\t{}, {} => {}", i, j, traceback[i][j].score);
        ops.push(traceback[i][j].op.clone());
        match traceback[i][j].op {
            Op::M => { },
            Op::D => { },
            Op::I => { },
        }
    }

    ops.reverse();

    traceback[m][n].score
}

