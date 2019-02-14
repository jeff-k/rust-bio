// Copyright 2019 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Coloured de Bruijn graph assembly into partially ordered graph

use alignment::poa;
use boomphf::hashmap::BoomHashMap2;
use debruijn::compression::*;
use debruijn::dna_string::{DnaString, DnaStringSlice};
use debruijn::filter::{filter_kmers, EqClassIdType};
use debruijn::graph::{BaseGraph, DebruijnGraph};
use debruijn::{kmer, Exts};
use petgraph::{Directed, Graph, Incoming};
use utils::TextSlice;

use std::str;

struct Debruijn {
    graph: BaseGraph<kmer::Kmer32, EqClassIdType>,
}

// copy+pasted from rust-pseudoaligner
lazy_static! {
    static ref PERM: Vec<usize> = {
        let maxp = 1 << (2 * 6);
        let mut permutation = Vec::with_capacity(maxp);
        for i in 0..maxp {
            permutation.push(i);
        }
        permutation
    };
}

impl Debruijn {
    pub fn new(tseq: TextSlice) -> Self {
        // the summarizer defines the criteria for admitting a kmer into the graph
        let summarizer = Box::new(debruijn::filter::CountFilterEqClass::new(1));
        let seq = DnaString::from_dna_string(str::from_utf8(tseq).unwrap());
        // ...
        let msps = debruijn::msp::simple_scan::<_, kmer::Kmer6>(32, &seq, &PERM, true);

        // ..
        let mut exts = vec![];
        for msp in msps {
            exts.push(Exts::from_dna_string(&seq, msp.start(), msp.end()));
        }

        // apply the summarizer to the hashmap
        let v = vec![(seq, exts, 0)];
        let (phf, _): (BoomHashMap2<kmer::Kmer32, Exts, EqClassIdType>, _) =
            filter_kmers(&v, &summarizer, true, true, 8);
        // filter_kmers returns (hashmap, Vec(kmer)) tuple

        Debruijn {
            // compress?
            graph: compress_kmers_with_hash(false, ScmapCompress::new(), &phf),
        }
    }

    pub fn push(&self, seq: TextSlice, seqid: u8) {}

    pub fn len(&self) -> usize {
        self.graph.len()
    }

    pub fn to_poa(
        &self,
        start: TextSlice,
        end: TextSlice,
    ) -> poa::Poa<poa::MatchFunc, i32, poa::EdgeFunc<i32>> {
        let mut graph: Graph<u8, i32, Directed, usize> = Graph::with_capacity(10000, 9999);

        // traverse the debruijn graph to produce the assembly
        //filter_kmers
        //       for node in self.graph.iter() {
        //            graph.add_node(petgraph::Node::new(node.to_ascii()));
        //       }
        let scoring = poa::Scoring::new(-1, 0, |a: u8, b: u8| {
            if a == b {
                1i32
            }
        });
        let labeller = poa::EdgeUpdater::new(1, |a: i32, b: i32| a + b);
        poa::Poa::new(scoring, labeller, graph)
    }

    pub fn to_gfa(&self) {
        //self.graph.to_gfa();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init_graph() {
        let g = Debruijn::new(b"ASDF");
        assert_eq!(g.len(), 0);
    }
}
