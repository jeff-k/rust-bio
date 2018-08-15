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

pub struct translation_table {
    translate: TextSlice -> TextSlice
}

impl translation_table {
    pub fn new(aas: TextSlice, starts: TextSlice, base1: TextSlice,
               base2: TextSlice, base3: TextSlice) -> translation_table {

        for (a, (b, c)) in base1.zip(base2.zip(base3)) {
            table[(a,b,c)] = protein;
        }
    }
}

pub fn translate(TextSlice) -> TextSlice {
}

// constructor for the standard genetic code from condensed translation table
// (NCBI Taxonomy group)

let table1 = translation_table::new(
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    b"---M---------------M---------------M----------------------------"
    b"TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    b"TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    b"TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG")

#[cfg(test)]
mod tests {
    use super::*;

}
