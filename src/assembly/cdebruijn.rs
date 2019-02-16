// Copyright 2019 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Coloured de Bruijn graph assembly into partially ordered graph

use alignment::pairwise::Scoring;
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
    pub fn new(tseq: TextSlice, useq: TextSlice) -> Self {
        // the summarizer defines the criteria for admitting a kmer into the graph
        let summarizer = Box::new(debruijn::filter::CountFilterEqClass::new(1));
        // ...

        let seqs = [
            DnaString::from_dna_string(str::from_utf8(tseq).unwrap()),
            DnaString::from_dna_string(str::from_utf8(useq).unwrap()),
        ];
        let msps = debruijn::msp::simple_scan::<_, kmer::Kmer6>(32, &seqs[0], &PERM, true);

        // ..
        //        let mut ext; // = vec![];
        //        let mut v = [];
        let mut phf: BoomHashMap2<kmer::Kmer32, Exts, EqClassIdType> = Default::default();
        //let msps = seq.iter().enumerate().map(|i, seq| { debruijn::msp::simple_scan::<_, kmer>(32, &seq, &PERM, true) }).collect();
        for msp in msps {
            //            ext = Exts::from_dna_string(&seq, msp.start(), msp.end());

            // i have no idea what im doing
            let seqq = DnaString::from_dna_string(str::from_utf8(tseq).unwrap());
            let seqqq = DnaString::from_dna_string(str::from_utf8(tseq).unwrap());
            let v = [(
                seqq,
                Exts::from_dna_string(&seqqq, msp.start(), msp.end()),
                0,
            )];
            // apply the summarizer to the hashmap
            let x = filter_kmers(&v, &summarizer, true, true, 8);
            phf = x.0;
            //return Debruijn {
            // compress?

            // filter_kmers returns (hashmap, Vec(kmer)) tuple
        }

        Debruijn {
            graph: compress_kmers_with_hash(false, ScmapCompress::new(), &phf),
        }
    }

    pub fn gfa(self, file: String) {
        self.graph.finish().to_gfa(file);
    }

    pub fn push(&self, seq: TextSlice, seqid: u8) {}

    pub fn len(&self) -> usize {
        self.graph.len()
    }

    pub fn to_poa(&self, start: TextSlice, end: TextSlice) {
        // -> poa::Poa<poa::MatchFunc, i32, poa::EdgeFunc<i32>> {
        let mut graph: Graph<u8, i32, Directed, usize> = Graph::with_capacity(10000, 9999);

        // traverse the debruijn graph to produce the assembly
        //filter_kmers
        //       for node in self.graph.iter() {
        //            graph.add_node(petgraph::Node::new(node.to_ascii()));
        //       }
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let labeller = poa::EdgeUpdater::new(1, |a: i32, b: i32| a + b);
        let p = poa::Poa::new(scoring, labeller, graph);
        p.write_dot("/tmp/poa_out.dot");
    }

    pub fn to_gfa(&self) {
        //self.graph.to_gfa();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::prelude::*;

    #[test]
    fn test_init_graph() {
        let mut fds = File::open("/tmp/hxb2_flat").unwrap();
        let mut contents = String::new();
        fds.read_to_string(&mut contents).unwrap();
        let g = Debruijn::go(contents.as_bytes());
        assert_eq!(0, 0);
    }
}
