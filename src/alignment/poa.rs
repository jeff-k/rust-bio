// Copyright 2017-2018 Brett Bowman, Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of multiple homologous sequences.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a modern reference implementation, see poapy:
//! https://github.com/ljdursi/poapy
//!
//! # Example
//!
//! ```
//! use bio::alignment::poa::*;
//! use bio::alignment::pairwise::Scoring;
//!
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
//! let mut aligner = Aligner::new(scoring, x, 1i32);
//! // z differs from x in 3 locations
//! assert_eq!(aligner.global(z).alignment().score, 1);
//! aligner.global(y).add_to_graph(1i32);
//! // z differs from x and y's partial order alignment by 1 base
//! assert_eq!(aligner.global(z).alignment().score, 5);
//! ```
//!

use std::cmp::{max, Ordering};
use std::fmt::Display;

use crate::utils::TextSlice;

use crate::alignment::pairwise::{MatchFunc, Scoring, Semiring};
use crate::alignment::AlignmentMode;

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

//pub type

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment<S: Semiring> {
    pub score: S,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
    mode: AlignmentMode,
}

#[derive(Debug, Clone)]
pub struct TracebackCell<S: Semiring> {
    score: S,
    op: AlignmentOperation,
}

//impl Default for TracebackCell { }

//impl Eq for TracebackCell {}

fn argmax2<S: Semiring>(t1: TracebackCell<S>, t2: TracebackCell<S>) -> TracebackCell<S> {
    if t1.score > t2.score {
        TracebackCell {
            score: t1.score.add(t2.score),
            op: t1.op,
        }
    } else if t2.score > t1.score {
        TracebackCell {
            score: t2.score.add(t1.score),
            op: t2.op,
        }
    } else {
        TracebackCell {
            score: t1.score,
            op: t2.op,
        }
    }
}

fn argmax3<S: Semiring>(
    t1: TracebackCell<S>,
    t2: TracebackCell<S>,
    t3: TracebackCell<S>,
) -> TracebackCell<S> {
    argmax2(t1, argmax2(t2, t3))
}

pub struct Traceback<S: Semiring + Copy + Clone + Display> {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<Vec<TracebackCell<S>>>,
}

impl<S: Semiring + Copy + Clone + Display> Traceback<S> {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    ///
    fn with_capacity(m: usize, n: usize) -> Self {
        let matrix = vec![
            vec![
                TracebackCell::<S> {
                    score: S::one(),
                    op: AlignmentOperation::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }

    /// Populate the edges of the traceback matrix
    fn initialize_scores(&mut self, gap_open: S) {
        for (i, row) in self
            .matrix
            .iter_mut()
            .enumerate()
            .take(self.rows + 1)
            .skip(1)
        {
            // TODO: pass a callback to customize clipping penalties
            row[0] = TracebackCell {
                score: S::one().mul(gap_open), // gap_open penalty
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=self.cols {
            self.matrix[0][j] = TracebackCell {
                score: S::one().mul(gap_open),
                op: AlignmentOperation::Ins(None),
            };
        }
    }

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell<S>) {
        self.matrix[i][j] = cell;
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell<S> {
        &self.matrix[i][j]
    }

    pub fn print(&self, g: &Graph<u8, S, Directed, usize>, query: TextSlice) {
        let (m, n) = (g.node_count(), query.len());
        print!(".\t");
        for base in query.iter().take(n) {
            print!("{:?}\t", *base);
        }
        for i in 0..m {
            print!("\n{:?}\t", g.raw_nodes()[i].weight);
            for j in 0..n {
                print!("{}.\t", self.get(i + 1, j + 1).score);
            }
        }
        println!();
    }

    pub fn alignment(&self) -> Alignment<S> {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.matrix[i][j].op.clone());
            match self.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Match(None) => {
                    break;
                }
                AlignmentOperation::Del(None) => {
                    j -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    i -= 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.matrix[self.last.index() + 1][self.cols].score,
            operations: ops,
            mode: AlignmentMode::Custom,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct Aligner<S: Semiring + Copy + Clone + Display, F: MatchFunc<S>> {
    sequence_names: Vec<String>,
    pub traceback: Traceback<S>,
    query: Vec<u8>,
    poa: Poa<S, F>,
}

impl<S: Semiring + Copy + Clone + Display, F: MatchFunc<S>> Aligner<S, F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<S, F>, reference: TextSlice, colour: S) -> Self {
        Aligner {
            sequence_names: vec![],
            traceback: Traceback::<S>::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(scoring, reference, colour),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self, colour: S) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, &self.query, colour);
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment<S> {
        self.traceback.alignment()
    }

    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: TextSlice) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(query);
        self
    }

    /// Return alignment graph.
    pub fn graph(&self) -> &Graph<u8, S, Directed, usize> {
        &self.poa.graph
    }
}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
///
pub struct Poa<S: Semiring, F: MatchFunc<S>> {
    scoring: Scoring<S, F>,
    pub graph: Graph<u8, S, Directed, usize>,
}

impl<S: Semiring + Copy + Clone + Display, F: MatchFunc<S>> Poa<S, F> {
    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `poa` - the partially ordered reference alignment
    ///
    pub fn new(scoring: Scoring<S, F>, graph: Graph<u8, S, Directed, usize>) -> Self {
        Poa { scoring, graph }
    }

    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    ///
    pub fn from_string(scoring: Scoring<S, F>, seq: TextSlice, colour: S) -> Self {
        let mut graph: Graph<u8, S, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, colour);
            prev = node;
        }

        Poa { scoring, graph }
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///
    pub fn global(&self, query: TextSlice) -> Traceback<S> {
        assert!(self.graph.node_count() != 0);

        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open);

        traceback.set(
            0,
            0,
            TracebackCell {
                score: S::one(),
                op: AlignmentOperation::Match(None),
            },
        );

        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    TracebackCell {
                        score: traceback
                            .get(0, j - 1)
                            .score
                            .mul(self.scoring.match_fn.score(r, *q)),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: S::zero(),
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        max_cell = argmax3(
                            max_cell,
                            TracebackCell {
                                score: traceback
                                    .get(i_p, j - 1)
                                    .score
                                    .mul(self.scoring.match_fn.score(r, *q)),
                                op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                            },
                            TracebackCell {
                                score: traceback.get(i_p, j).score.mul(self.scoring.gap_open),
                                op: AlignmentOperation::Del(Some((i_p - 1, i))),
                            },
                        );
                    }
                    max_cell
                };

                let score = argmax2(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score.mul(self.scoring.gap_open),
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
            }
        }

        traceback
    }

    /// Experimental: return sequence of traversed edges
    ///
    /// Only supports alignments for sequences that have already been added,
    /// so all operations must be Match.
    pub fn edges(&self, aln: Alignment<S>) -> Vec<usize> {
        let mut path: Vec<usize> = vec![];
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut _i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    _i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    let edge = self.graph.find_edge(prev, node).unwrap();
                    path.push(edge.index());
                    prev = NodeIndex::new(p);
                    _i += 1;
                }
                AlignmentOperation::Ins(None) => {}
                AlignmentOperation::Ins(Some(_)) => {}
                AlignmentOperation::Del(_) => {}
            }
        }
        path
    }

    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    ///
    pub fn add_alignment(&mut self, aln: &Alignment<S>, seq: TextSlice, colour: S) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    if seq[i] != self.graph.raw_nodes()[*p].weight {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, colour);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                let n_c = self.graph.edge_weight(edge).unwrap().mul(colour);
                                self.graph.update_edge(prev, node, n_c);
                            }
                            None => {
                                // where the previous node was newly added
                                self.graph.add_edge(prev, node, colour);
                            }
                        }
                        prev = NodeIndex::new(*p);
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, colour);
                    prev = node;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::pairwise::Scoring;
    use petgraph::graph::NodeIndex;

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph

        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let poa = Poa::from_string(scoring, b"123456789", 1i32);
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let poa = Poa::from_string(scoring, b"GATTACA", 1i32);
        let alignment = poa.global(b"GCATGCU").alignment();
        assert_eq!(alignment.score, 0);

        let alignment = poa.global(b"GCATGCUx").alignment();
        assert_eq!(alignment.score, -1);

        let alignment = poa.global(b"xCATGCU").alignment();
        assert_eq!(alignment.score, -2);
    }

    #[test]
    fn test_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let seq1 = b"TTTTT";
        let seq2 = b"TTATT";
        let mut poa = Poa::from_string(scoring, seq1, 1i32);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2).alignment();
        assert_eq!(alignment.score, 3);
    }

    #[test]
    fn test_alt_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCTTAA";
        let seq2 = b"TTTTGGAA";
        let mut poa = Poa::from_string(scoring, seq1, 1i32);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2).alignment();
        poa.add_alignment(&alignment, seq2, 1i32);
        assert_eq!(poa.graph.edge_count(), 14);
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(5), NodeIndex::new(10)));
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(11), NodeIndex::new(6)));
    }

    #[test]
    fn test_insertion_on_branch() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTATGGGAA";
        let seq3 = b"TTGGTTTGCGAA";
        let mut poa = Poa::from_string(scoring, seq1, 1i32);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        let node3 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, node3, 1);
        poa.graph.add_edge(node3, tail, 1);
        let alignment = poa.global(seq2).alignment();
        assert_eq!(alignment.score, 2);
        poa.add_alignment(&alignment, seq2, 1i32);
        let alignment2 = poa.global(seq3).alignment();

        assert_eq!(alignment2.score, 10);
    }

    #[test]
    fn test_poa_method_chaining() {
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"TTCCGGTTTAA", 1i32);
        aligner
            .global(b"TTGGTATGGGAA")
            .add_to_graph(1i32)
            .global(b"TTGGTTTGCGAA")
            .add_to_graph(1i32);
        assert_eq!(aligner.alignment().score, 10);
    }
}
