// Copyright 2019 Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Coloured de Bruijn graph assembly

use debruijn::filter::EqClassIdType;
use debruijn::graph::BaseGraph;
use debruijn::kmer;
use utils::TextSlice;

struct Debruijn {
    graph: BaseGraph<kmer::Kmer24, EqClassIdType>,
}

impl Debruijn {
    pub fn new(seq: TextSlice) -> Self {
        Debruijn {
            graph: BaseGraph::new(false),
        }
    }

    pub fn len(&self) -> usize {
        self.graph.len()
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
