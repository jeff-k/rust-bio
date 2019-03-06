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
use debruijn::{kmer, Dir, Exts, Kmer};
use petgraph::{Directed, Graph, Incoming};
use utils::TextSlice;

use std::str;

struct Debruijn {
    graph: DebruijnGraph<kmer::Kmer32, EqClassIdType>,
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
        let msps = debruijn::msp::simple_scan::<_, kmer::Kmer6>(32, &seq, &PERM, true);

        let mut v: Vec<_> = vec![];
        for msp in msps {
            v.push((
                DnaString::from_dna_string(str::from_utf8(tseq).unwrap()),
                Exts::from_dna_string(&seq, msp.start(), msp.end()),
                0,
            ));
        }
        // apply the summarizer to the hashmap
        let x = filter_kmers(&v, &summarizer, true, true, 8);
        let phf: BoomHashMap2<kmer::Kmer32, Exts, EqClassIdType> = x.0;

        Debruijn {
            graph: compress_kmers_with_hash(true, ScmapCompress::new(), &phf).finish(),
        }
    }

    pub fn to_gfa(self, file: String) {
        self.graph.to_gfa(file);
    }

    pub fn push(&self, seq: TextSlice, seqid: u8) {}

    pub fn len(&self) -> usize {
        self.graph.len()
    }

    pub fn to_poa(&self, start: TextSlice, end: TextSlice) {
        //-> Option(poa::Poa<poa::MatchFunc, i32, poa::EdgeFunc<i32>>) {
        //let mut graph: Graph<u8, i32, Directed, usize> = Graph::new();

        // build a partial order from a list of paths that begin and end in the same place
        let snode = match self
            .graph
            .search_kmer(kmer::Kmer32::from_ascii(start), Dir::Left)
        {
            Some(idx) => idx,
            None => 0,
        };

        let enode = match self
            .graph
            .search_kmer(kmer::Kmer32::from_ascii(end), Dir::Right)
        {
            Some(idx) => idx,
            None => 0,
        };

        let path = self.graph.max_path(|_| 1.0, |_| true);
        let seq = self.graph.sequence_of_path(path.iter());
        println!("node id {:?}/{:?} for end kmer: {:?}", enode, self.graph.len(), kmer::Kmer32::from_ascii(end).to_string());
        println!("node id {:?}/{:?} for start kmer: {:?}", snode, self.graph.len(), kmer::Kmer32::from_ascii(start).to_string());
        let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
        let labeller = poa::EdgeUpdater::new(1, |a: i32, b: i32| a + b);
        let p: poa::Poa<_, i32, _> = poa::Poa::from_string(scoring, labeller, &seq.to_ascii_vec());

        p.write_dot("/tmp/poa_out.dot");

        //poa
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
        let g = Debruijn::new(contents.as_bytes());
        g.to_gfa("/tmp/graph_out".to_string());
        assert_eq!(0, 0);
    }
//    #[test]
    fn test_to_poa_again() {
        let mut fds = File::open("/tmp/hxb2_flat").unwrap();
        let mut contents = String::new();
        fds.read_to_string(&mut contents).unwrap();
        let g = Debruijn::new(contents.as_bytes());
        g.to_poa(b"TGGAAGGGCTAATTCACTCCCAACGAAGACAA",
                 b"ATCAGATATCCACTGACCTTTGGATGGTGCTA");
        assert_eq!(0, 0);
    }
//    #[test]
    fn test_to_poa_insidemer() {
        let mut fds = File::open("/tmp/hxb2_flat").unwrap();
        let mut contents = String::new();
        fds.read_to_string(&mut contents).unwrap();
        let g = Debruijn::new(contents.as_bytes());
        g.to_poa(b"GAAGGGCTAATTCACTCCCAACGAAGACAAGA",
                 b"CTGATTAGCAGAACTACACACCAGGGCCAGGG");
        assert_eq!(0, 0);
    }


//    #[test]
    fn test_to_poa() {
        let mut fds = File::open("/tmp/hxb2_flat").unwrap();
        let mut contents = String::new();
        fds.read_to_string(&mut contents).unwrap();
        let g = Debruijn::new(contents.as_bytes());
        g.to_poa(b"ACGAAGACAAGATATCCTTGATCTGTGGATCT",
                 b"ATGGCCCGAGAGCTGCATCCGGAGTACTTCAA");
        assert_eq!(0, 0);
    }
}
