// Copyright 2018 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Frameshift detection by partial order alignment
//!

use std::cmp::{max, Ordering};

use utils::TextSlice;

use alignment::pairwise::{MatchFunc, Scoring};
use alignment::AlignmentMode;

use alphabets::translation::{nuc2amino, table1, Translation_Table};

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
    Frameshift(Option<(usize, usize)>),
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
    fs_penalty: i32,
    table: Translation_Table,
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
    pub fn from_string(
        reference: TextSlice,
        table: Translation_Table,
        fs_penalty: i32,
        gap_open: i32,
        match_fn: F,
    ) -> Self {
        let s = vec![b'^', b'$'];
        let d = [&s[0..1], &nuc2amino(reference, &table), &s[1..2]].concat();
        let mut graph: Graph<u8, i32, Directed, usize> =
            Graph::with_capacity(d.len() * 3, d.len() * 3);
        let mut prev: NodeIndex<usize> = graph.add_node(d[0]);
        let mut node: NodeIndex<usize>;
        for i in 1..d.len() {
            node = graph.add_node(d[i]);
            graph.add_edge(prev, node, 0);
            prev = node;
        }

        let mut fs1p: NodeIndex<usize> = NodeIndex::new(0);
        let mut fs1pp: NodeIndex<usize> = NodeIndex::new(0);

        let mut fs2p: NodeIndex<usize> = NodeIndex::new(0);
        let mut fs2pp: NodeIndex<usize> = NodeIndex::new(0);

        let mut prev: NodeIndex<usize> = NodeIndex::new(0);

        //        let mut fs1_2: NodeIndex<usize> = NodeIndex::new(0);
        //        let mut fs2_2: NodeIndex<usize> = NodeIndex::new(0);

        let s1 = nuc2amino(&reference[1..], &table);
        let s2 = nuc2amino(&reference[2..], &table);
        let last = NodeIndex::new(graph.node_count() - 1);
        for i in 0..(graph.node_count() - 2) {
            let node = NodeIndex::new(i + 1);

            let fs1 = graph.add_node(s1[i]);
            let fs2 = graph.add_node(s2[i]);

            if i > 1 {
                graph.add_edge(fs1pp, node, 2 * fs_penalty);
                graph.add_edge(fs2pp, node, 1 * fs_penalty);
                graph.add_edge(fs2pp, fs1, 2 * fs_penalty);
            }

            if i > 0 {
                //graph.add_edge(fs1p, node, -2);
                graph.add_edge(fs1p, fs2, 1 * fs_penalty);
                graph.add_edge(fs1p, fs1, 0);
                graph.add_edge(fs2p, fs2, 0);
            }

            graph.add_edge(prev, fs1, 1 * fs_penalty);
            graph.add_edge(prev, fs2, 2 * fs_penalty);

            fs1pp = fs1p;
            fs2pp = fs2p;

            fs1p = fs1;
            fs2p = fs2;
            prev = node;
        }
        graph.add_edge(fs1p, last, 0);
        graph.add_edge(fs2p, last, 0);

        Aligner {
            scoring: Scoring::new(gap_open, 0, match_fn),
            fs_penalty: fs_penalty,
            table: table,
            traceback: Traceback::with_capacity(0, 0),
            graph: graph,
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
                        let edge = self.graph.find_edge(prevs[prev_n], node).unwrap();
                        let weight = self.graph.raw_edges()[edge.index()].weight;
                        if weight == 0 {
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
                        } else {
                            max_cell = max(
                                max_cell,
                                TracebackCell {
                                    score: self.traceback.matrix[i_p][j - 1].score + weight,
                                    op: AlignmentOperation::Frameshift(Some((i_p - 1, i))),
                                },
                            );
                        };
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
                AlignmentOperation::Frameshift(None) => {}
                AlignmentOperation::Frameshift(Some((p, _))) => {
                    i = p + 1;
                    j = j - 1;
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

    pub fn reconstruct(&self, seq: TextSlice, aln: Alignment) -> Vec<u8> {
        //        let mut aln = self.global(seq);
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        let mut out = Vec::new();

        for op in aln.operations.iter() {
            if i + 3 >= seq.len() {
                break;
            }
            let codon = &seq[i..i + 3];
            println!("{:?} - {:?}", String::from_utf8_lossy(codon), op);
            match op {
                AlignmentOperation::Match(None) => {
                    i = i + 3;
                }
                AlignmentOperation::Match(Some((_, _))) => {
                    out.push(codon);
                    i = i + 3;
                }
                AlignmentOperation::Frameshift(_) => {
                    out.push(b"nnn");
                    i = i + 1;
                }

                AlignmentOperation::Ins(_) => {
                    i = i + 3;
                    out.push(codon);
                }
                AlignmentOperation::Del(_) => {
                    i = i + 3;
                    out.push(codon);
                }
            }
        }
        out.concat()
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
    use alphabets::translation::{nuc2amino, table1, test_table, Translation_Table};
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
    fn test_braid() {
        let match_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut braid =
            Aligner::from_string(b"abcdefghijklmnopqrxx", test_table(), -1, -2, match_fn);
        let aln = braid.global(b"abcdefghiijklmnopqr");
        println!(
            "reconstruction {:?}",
            String::from_utf8_lossy(&braid.reconstruct(b"abcdefghiijklmnopqr", aln))
        );
        //assert!(false);
    }

}
