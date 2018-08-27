// Copyright 2017 Brett Bowman
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of long homologous sequences.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a reference implementation that inspired this code, see poapy:
//! https://github.com/ljdursi/poapy

use std::cmp::{max, Ordering};
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::str;
use std::string::String;

use alphabets::translation::table1;
use utils::TextSlice;

use petgraph::dot::Dot;
use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;
use petgraph::{Directed, Graph, Incoming};

pub const MIN_SCORE: i32 = -858_993_459; // negative infinity; see alignment/pairwise/mod.rs

/// A Partial-Order Alignment Graph
pub struct POAGraph {
    // a POAGraph struct stores a reference DAG and labels the edges with the
    // count or names of the input sequences (TODO)
    graph: Graph<u8, i32, Directed, usize>,
}

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
#[derive(Debug, Clone)]
pub enum Op {
    Match(Option<(usize, usize)>),
    //Del(Option<(usize, usize)>),
    //Ins(Option<usize>),
    Fs(Option<(usize, usize)>),
}

pub struct Alignment {
    score: i32,
    operations: Vec<Op>,
}

pub struct Aligner {
    scoring: fn(u8, u8) -> i32,
}

#[derive(Debug, Clone)]
struct TracebackCell {
    score: i32,
    op: Op,
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

impl Eq for TracebackCell {}

impl Aligner {
    pub fn new() -> Self {
        let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        Aligner { scoring: score_fn }
    }

    /// Naive Needleman-Wunsch
    /// Populates the traceback matrix in $O(n^2)$ time using
    /// petgraph's constant time topological traversal
    pub fn global(&mut self, g: &Graph<u8, i32, Directed, usize>, query: TextSlice) -> Alignment {
        // dimensions of the traceback matrix
        let (m, n) = (g.node_count(), query.len());

        // initialize matrix
        let mut traceback: Vec<Vec<TracebackCell>> = vec![
            vec![
                TracebackCell {
                    score: 0,
                    op: Op::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        let mut ops: Vec<Op> = vec![];

        for i in 1..(m + 1) {
            // TODO: these should be -1 * distance from head node
            traceback[i][0] = TracebackCell {
                score: -1 * i as i32,
                op: Op::Match(None),
            };
        }
        for j in 1..(n + 1) {
            traceback[0][j] = TracebackCell {
                score: -1 * j as i32,
                op: Op::Match(None),
            };
        }

        traceback[0][0] = TracebackCell {
            score: 0,
            op: Op::Match(None),
        };

        // construct the score matrix (naive)
        let mut topo = Topo::new(&g);

        // store the last visited node in topological order so that
        // we can index into the end of the alignment when we backtrack
        let mut last: NodeIndex<usize> = NodeIndex::new(0);
        while let Some(node) = topo.next(&g) {
            // reference base and index
            let r = g.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> = g.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            let mut j_p: usize = 0;
            for (j_p, q) in query.iter().enumerate() {
                //while j_p < query.len() - 1 {
                let j = j_p + 1;

                //                for fs in -1..1 {
                //                    let &protein =
                //                        table1().translate(&query[(j as i8 + fs) as usize..(j + fs as usize) + 3]);
                //                        //String::from_utf8_lossy(&[protein])
                //                    );
                //                    let op = if fs == 0 { Op::Match(None) } else { Op::Fs(fs) };
                //                }

                // match and frameshift scores for the first reference base
                let (mat, fs) = if prevs.len() == 0 {
                    (
                        TracebackCell {
                            score: traceback[0][j - 1].score + (self.scoring)(r, *q),
                            op: Op::Match(None),
                        },
                        TracebackCell {
                            score: traceback[0][j].score - 1i32,
                            op: Op::Fs(None),
                        },
                    )
                } else {
                    let mut mat_max = TracebackCell {
                        score: MIN_SCORE,
                        op: Op::Match(None),
                    };
                    let mut fs_max = TracebackCell {
                        score: MIN_SCORE,
                        op: Op::Fs(None),
                    };

                    // iterate over previous nodes and their edge weights
                    for prev_n in 0..prevs.len() {
                        let i_p: usize = prevs[prev_n].index() + 1; // index of previous node
                        let edge = g.find_edge(prevs[prev_n], node).unwrap();
                        let weight = g.raw_edges()[edge.index()].weight;
                        if weight == 0 {
                            mat_max = max(
                                mat_max,
                                TracebackCell {
                                    score: traceback[i_p][j - 1].score + (self.scoring)(r, *q),
                                    //                                op: Op::Match(Some((i_p - 1, i - 1))),
                                    op: Op::Match(Some((i_p - 1, i - 1))),
                                },
                            );
                        } else {
                            fs_max = max(
                                fs_max,
                                TracebackCell {
                                    score: traceback[i_p][j - 1].score + weight,
                                    op: Op::Fs(Some((i_p - 1, i))),
                                },
                            );
                        }
                    }
                    (mat_max, fs_max)
                };
                let score = max(mat, fs);
                traceback[i][j] = score;
            }
        }

        //        dump_traceback(&traceback, g, query);

        // Now backtrack through the matrix to construct an optimal path
        let mut i = last.index() + 1;
        let mut j = n;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(traceback[i][j].op.clone());
            match traceback[i][j].op {
                Op::Match(Some((p, _))) => {
                    i = p + 1;
                    j = j - 1;
                }
                Op::Match(None) => {
                    break;
                }
                Op::Fs(None) => {}
                Op::Fs(Some((p, _))) => {
                    i = p + 1;
                    j = j - 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: traceback[last.index() + 1][n].score,
            operations: ops,
        }
    }
}

impl POAGraph {
    /// Create new POAGraph instance initialized with a sequence
    ///
    /// # Arguments
    ///
    /// * `label` - sequence label for an initial sequence
    /// * `sequence` - TextSlice from which to initialize the POA
    ///
    pub fn new(_label: &str, sequence: TextSlice) -> POAGraph {
        let graph = POAGraph::seq_to_graph(sequence);
        POAGraph { graph: graph }
    }

    /// Add a new unaligned sequence to the underlying graph.
    /// Useful for both initializing the graph with its first sequence, as
    /// well as adding unaligned prefix or suffix sequence from partially
    /// aligned reads.
    ///
    /// # Arguments
    ///
    /// * `label` - the id of the sequence to be added
    /// * `sequence` - The sequence to be added
    ///
    fn seq_to_graph(sequence: TextSlice) -> Graph<u8, i32, Directed, usize> {
        // this should return a POAGraph
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(sequence.len(), sequence.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(sequence[0]);
        let mut node: NodeIndex<usize>;
        for i in 1..sequence.len() {
            node = graph.add_node(sequence[i]);
            graph.add_edge(prev, node, 0);
            prev = node;
        }
        graph
    }

    /// Align a sequence to the current graph and return the
    /// resulting alignment object
    ///
    /// # Arguments
    ///
    /// * `sequence` - The new sequence to align as a text slice
    ///
    pub fn align_sequence(&self, seq: &[u8]) -> Alignment {
        let mut aligner = Aligner::new();
        aligner.global(&self.graph, seq)
    }

    /// Incorporate a new sequence into the graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment object of the new sequence to the graph
    /// * `label` - The name of the new sequence being added to the graph
    /// * `seq` - The complete sequence of the read being incorporated
    ///
    pub fn incorporate_alignment(&mut self, aln: Alignment, _label: &str, seq: TextSlice) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations {
            match op {
                Op::Match(None) => {
                    i = i + 1;
                }
                Op::Match(Some((_, p))) => {
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
                Op::Fs(_) => {}
            }
        }
    }

    pub fn braid(&mut self, s1: TextSlice, s2: TextSlice) {
        let mut fs1p: NodeIndex<usize> = NodeIndex::new(0);
        let mut fs1pp: NodeIndex<usize> = NodeIndex::new(0);

        let mut fs2p: NodeIndex<usize> = NodeIndex::new(0);
        let mut fs2pp: NodeIndex<usize> = NodeIndex::new(0);

        let mut prev: NodeIndex<usize> = NodeIndex::new(0);

        //        let mut fs1_2: NodeIndex<usize> = NodeIndex::new(0);
        //        let mut fs2_2: NodeIndex<usize> = NodeIndex::new(0);

        let last = NodeIndex::new(self.graph.node_count() - 1);
        for i in 0..(self.graph.node_count() - 2) {
            let node = NodeIndex::new(i + 1);

            let fs1 = self.graph.add_node(s1[i]);
            let fs2 = self.graph.add_node(s2[i]);

            if i > 1 {
                self.graph.add_edge(fs2pp, node, -2);
                self.graph.add_edge(fs2pp, fs1, -2);
            }

            if i > 0 {
                self.graph.add_edge(fs1p, node, -1);
                self.graph.add_edge(fs1p, fs2, -1);
                self.graph.add_edge(fs1p, fs1, 0);
                self.graph.add_edge(fs2p, fs2, 0);
            }

            self.graph.add_edge(prev, fs1, -1);
            self.graph.add_edge(prev, fs2, -2);

            fs1pp = fs1p;
            fs2pp = fs2p;

            fs1p = fs1;
            fs2p = fs2;
            prev = node;
        }
        self.graph.add_edge(fs1p, last, 0);
        self.graph.add_edge(fs2p, last, 0);
    }
    /// Write the current graph to a specified filepath in dot format for
    /// visualization, primarily for debugging / diagnostics
    ///
    /// # Arguments
    ///
    /// * `filename` - The filepath to write the dot file to, as a String
    ///
    pub fn write_dot(&self, filename: String) {
        let mut file = match File::create(&filename) {
            Err(why) => panic!("couldn't open file {}: {}", filename, why.description()),
            Ok(file) => file,
        };
        let g = self.graph.map(|_, nw| *nw as char, |_, ew| ew);
        match file.write_all(Dot::new(&g).to_string().as_bytes()) {
            Err(why) => panic!("couldn't write to file {}: {}", filename, why.description()),
            _ => (),
        }
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
        print!("{:?}\t", String::from_utf8_lossy(&[query[i]]));
    }
    for i in 0..m {
        print!(
            "\n{:?}\t",
            String::from_utf8_lossy(&[g.raw_nodes()[i].weight])
        );
        for j in 0..n {
            print!("{}.\t", traceback[i + 1][j + 1].score);
        }
    }
    print!("\n");
}

#[cfg(test)]
mod tests {
    use alignment::poa::POAGraph;
    use alphabets::translation::table1;
    use petgraph::graph::NodeIndex;

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
    fn test_frameshift() {
        //        let dna1 = b"ACGTGCGGAATCGCGACGTGGGTAGC";
        let dna1 = b"^TCGIATWV$";
        let dna2 = b"RAESRRG*";
        let dna3 = b"VRNRDVGS";
        //        let f1 = table1().translate(dna1);
        //        let f2 = table1().translate(dna2);
        //        let f3 = table1().translate(dna3);

        let mut poa = POAGraph::new("seq", dna1);
        poa.braid(dna2, dna3);
        //        poa.write_dot("/tmp/x.dot".to_string());

        let test = b"^TCGIGRG*$";
        let alignment = poa.align_sequence(test);
        //        println!("{:?}", alignment.operations);
        //        println!("final score: {:?}", alignment.score);
        //        assert!(false);
    }

}
