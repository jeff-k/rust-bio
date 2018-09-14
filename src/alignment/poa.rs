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
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
//! let mut aligner = Aligner::from_string(x, -1, &score);
//! assert_eq!(aligner.global(z).score, 1);
//! aligner.add_sequence(y);
//! // z differs from its closest homologue in the graph by one base
//! assert_eq!(aligner.global(z).score, 5);
//! ```

use std::cmp::{max, Ordering};

use utils::TextSlice;

use alignment::pairwise::{MatchFunc, Scoring};
use alignment::AlignmentMode;

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
    Fs(Option<(usize, i32)>),
}

pub struct Alignment {
    pub score: i32,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
    mode: AlignmentMode,
}

#[derive(Debug, Clone)]
struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<Vec<TracebackCell>>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    ///
    fn with_capacity(m: usize, n: usize) -> Self {
        let mut matrix = vec![
            vec![
                TracebackCell {
                    score: 0,
                    op: AlignmentOperation::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        for i in 1..(m + 1) {
            // TODO: these should be -1 * distance from head node
            matrix[i][0] = TracebackCell {
                score: -1 * i as i32, // gap_open penalty
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..(n + 1) {
            matrix[0][j] = TracebackCell {
                score: -1 * j as i32,
                op: AlignmentOperation::Ins(None),
            };
        }

        Traceback {
            rows: m,
            cols: n,
            matrix: matrix,
        }
    }
}

/// A global aligner on partially ordered graphs
///
/// Internally stores a directed acyclic graph datastructure that informs the topology of the
/// traceback matrix. A mutable Aligner may have additional sequences added to the internal
/// partially ordered graph.
///
/// The partial order is expressed with a directed acyclic graph and informs the traversal of the
/// traceback matrix.
///
pub struct Aligner<F: MatchFunc> {
    scoring: Scoring<F>,
    traceback: Traceback,
    pub graph: Graph<u8, i32, Directed, usize>,
}

impl<F: MatchFunc> Aligner<F> {
    /// Create a new aligner instance from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    /// * `gap_open` - the negative score assigned when branching from the reference graph
    /// * `match_fn` - the pairwise score for substitutions (see bio::scores)
    ///
    pub fn from_string(reference: TextSlice, gap_open: i32, match_fn: F) -> Self {
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(reference.len(), reference.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(reference[0]);
        let mut node: NodeIndex<usize>;
        for i in 1..reference.len() {
            node = graph.add_node(reference[i]);
            graph.add_edge(prev, node, 1);
            prev = node;
        }

        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            traceback: Traceback::with_capacity(0, 0),
            graph: graph,
        }
    }

    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `poa` - the partially ordered reference alignment
    /// * `gap_open` - the negative score assigned when branch from the reference graph
    /// * `match_fn` - the pairwise score for substitutions (see bio::scores)
    ///
    pub fn from_graph(poa: Graph<u8, i32, Directed, usize>, gap_open: i32, match_fn: F) -> Self {
        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            traceback: Traceback::with_capacity(0, 0),
            graph: poa,
        }
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///
    pub fn global(&mut self, query: TextSlice) -> Alignment {
        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        self.traceback = Traceback::with_capacity(m, n);

        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        self.traceback.matrix[0][0] = TracebackCell {
            score: 0,
            op: AlignmentOperation::Match(None),
        };

        // construct the score matrix (naive)
        let mut topo = Topo::new(&self.graph);

        // store the last visited node in topological order so that
        // we can index into the end of the alignment when we backtrack
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.len() == 0 {
                    TracebackCell {
                        score: self.traceback.matrix[0][j - 1].score
                            + self.scoring.match_fn.score(r, *q),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: self.traceback.matrix[i_p][j - 1].score
                                        + self.scoring.match_fn.score(r, *q),
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: self.traceback.matrix[i_p][j].score
                                        + self.scoring.gap_open,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let score = max(
                    max_cell,
                    TracebackCell {
                        score: self.traceback.matrix[i][j - 1].score + self.scoring.gap_open,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                self.traceback.matrix[i][j] = score;
            }
        }

        //dump_traceback(&self.traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.traceback.matrix[i][j].op.clone());
            match self.traceback.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j = j - 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j = j - 1;
                }
                AlignmentOperation::Match(None) => {
                    break;
                }
                AlignmentOperation::Del(None) => {
                    j = j - 1;
                }
                AlignmentOperation::Ins(None) => {
                    i = i - 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.traceback.matrix[last.index() + 1][n].score,
            operations: ops,
            mode: AlignmentMode::Custom,
        }
    }

    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    ///
    pub fn align_sequence(&self, seq: &[u8]) -> Alignment {
        let mut aligner = Aligner::new();
        aligner.global(&self.graph, seq)
    }

    pub fn reconstruct(&self, aln: Alignment, _label: &str, seq: TextSlice) -> Vec<u8> {
=======
    pub fn add_alignment(&mut self, aln: Alignment, seq: TextSlice) {
>>>>>>> partial-order-aligner
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    i = i + 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    if seq[i] != self.graph.raw_nodes()[p].weight {
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                // where the previous node was newly added
                                self.graph.add_edge(prev, node, 1);
                            }
                        }
                        prev = NodeIndex::new(p);
                    }
                    i = i + 1;
                }
                AlignmentOperation::Ins(None) => {
                    i = i + 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    self.graph.add_edge(prev, node, 1);
                    prev = node;
                    i = i + 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }

    /// Align and incorporate a sequence to the partially ordered graph from a TextSlice
    ///
    /// # Arguments
    /// * `seq` - TextSlice to incorporate into the partial order aligner's graph
    ///
    pub fn add_sequence(&mut self, seq: TextSlice) {
        let alignment = self.global(seq);
        self.add_alignment(alignment, seq);
    }
}

// print out a traceback matrix
#[allow(dead_code)]
fn dump_traceback(
    traceback: &Vec<Vec<TracebackCell>>,
    g: &Graph<u8, i32, Directed, usize>,
    query: TextSlice,
) {
    let (m, n) = (g.node_count(), query.len());
    print!(".\t");
    for i in 0..n {
        print!("{:?}\t", query[i]);
    }
    for i in 0..m {
        print!("\n{:?}\t", g.raw_nodes()[i].weight);
        for j in 0..n {
            print!("{}.\t", traceback[i + 1][j + 1].score);
        }
    }
    print!("\n");
}

#[cfg(test)]
mod tests {
    use alignment::poa::Aligner;
    use petgraph::dot::Dot;
    use petgraph::graph::NodeIndex;
    use petgraph::{Directed, Graph};
    use std::error::Error;
    use std::fs::File;
    use std::io::Write;

    /// Write the current graph to a specified filepath in dot format for
    /// visualization, primarily for debugging / diagnostics
    ///
    /// # Arguments
    ///
    /// * `filename` - The filepath to write the dot file to, as a String
    ///
    #[allow(dead_code)]
    fn write_dot(graph: Graph<u8, i32, Directed, usize>, filename: String) {
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't open file {}: {}", filename, why.description()),
            Ok(file) => file,
        };
        let g = graph.map(|_, nw| *nw as char, |_, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
    }

    #[test]
    fn test_braiding() {
        let dna1 = b"^ABC";
        let dna2 = b"FGH";
        let dna3 = b"KLM";

        let mut poa = POAGraph::new("seq", dna1);
        poa.braid(dna2, dna3);
        //        poa.write_dot("/tmp/y.dot".to_string());
    }

    #[test]
    fn construct_frames() {
        let seq = b"abcdefghijklmnopqrxx";
        println!(
            "translated {:?}",
            String::from_utf8_lossy(&nuc2amino(seq, &test_table()))
        );
        let poa = braid_fs(seq, &test_table());
        let tst = b"abcdefghiijklmnopqr";
        println!(
            "correcting {:?}",
            String::from_utf8_lossy(&nuc2amino(tst, &test_table()))
        );
        let alignment = poa.align_sequence(&nuc2amino(tst, &test_table()));
        poa.write_dot("/tmp/test_braid.dot".to_string());
        println!("score {:?}", alignment.score);
        let r = poa.reconstruct(alignment, "asdf", tst);
        println!("{:?}", String::from_utf8_lossy(&r));
        println!(
            "{:?}",
            String::from_utf8_lossy(&nuc2amino(&r, &test_table()))
        );
        assert_eq!(r, b"abcdefghijklmnopqr");
    }

    #[test]
    fn test_simple_fs() {
        // ref
        let seq = b"ACGTGCGGATCGCGA";

        // query
        let dna = b"ACGTGCGGATCGCGA";

        let poa = braid_fs(dna, &table1());
        poa.write_dot("/tmp/simple.dot".to_string());
        //        let seq = b"ACGTGCATCGCGANN";
        let tst = &nuc2amino(seq, &table1());
        let s = vec![b'^', b'$'];
        let x = &[&s[0..1], tst, &s[1..2]].concat();
        println!("frameshifting {:?}", String::from_utf8_lossy(x));
        let aln = poa.align_sequence(x);
        println!("scored {:?}", aln.score);
        let r = poa.reconstruct(aln, "asdf", dna);
        println!("reconstruction: {:?}", String::from_utf8_lossy(&r));
        println!("{:?}", String::from_utf8_lossy(&nuc2amino(&r, &table1())));
        assert!(false);
    }
}
