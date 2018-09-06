// Copyright 2018 jeff-k
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Translate DNA to amino acid sequences
//!
//! # Example
//!
//! ```
//! ```

use std::collections::HashMap;
use utils::TextSlice;

pub struct Translation_Table {
    //    translate: TextSlice -> TextSlice,
    table: HashMap<(u8, u8, u8), u8>,
}

impl Translation_Table {
    pub fn new(
        aas: TextSlice,
        starts: TextSlice,
        base1: TextSlice,
        base2: TextSlice,
        base3: TextSlice,
    ) -> Translation_Table {
        let mut table = HashMap::new();
        for i in 0..aas.len() {
            table.insert((base1[i], base2[i], base3[i]), aas[i]);
        }
        Translation_Table { table: table }
    }

    pub fn translate(&self, ns: &[u8]) -> &u8 {
        match self.table.get(&(ns[0], ns[1], ns[2])) {
            Some(x) => x,
            None => &b'X',
        }
    }
}

pub fn nuc2amino(dna: TextSlice, table: &Translation_Table) -> Vec<u8> {
    // pad sequence with Ns
    let mut fs: Vec<u8> = vec![b'X'; (dna.len() / 3) + 1];
    for (i, codon) in dna.chunks(3).enumerate() {
        if codon.len() != 3 {
            fs[i] = b'X';
        } else {
            fs[i] = *table.translate(codon);
        }
    }
    fs
}

pub fn table1() -> Translation_Table {
    Translation_Table::new(
        b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        b"---M---------------M---------------M----------------------------",
        b"TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
        b"TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
        b"TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG",
    )
}

pub fn test_table() -> Translation_Table {
    Translation_Table::new(
        b"ABCABCDEFGHIJKLMNOPQ",
        b"---------",
        b"abcabcdefghijklmnopq",
        b"abcbcdefghijklmnopqr",
        b"abccdefghijklmnopqrs",
        )
}

// constructor for the standard genetic code from condensed translation table
// (NCBI Taxonomy group)

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translation() {
        let t = table1();
        assert_eq!(t.translate(b"TTT"), &b'F');
        assert_ne!(t.translate(b"GGT"), &b'T');
        assert_eq!(t.translate(b"GGT"), &b'G');
    }
}
