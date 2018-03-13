// Template Copyright 2014-2017 M. Rizky Luthfianto.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

// Empirical score matrix based on 25% divergent HIV sequences
// Nickle, David C., et al.
// "HIV-specific probabilistic models of protein evolution."
// PLoS One 2.6 (2007): e503.

use ndarray;


lazy_static! {
    // taken from https://github.com/cfe-lab/MiCall/blob/master/micall/alignment/gotoh.cpp
	static ref MAT: ndarray::Array2<i32> = ndarray::Array::from_shape_vec((24, 24), vec![
        7,-7,-7,-4,-10,-11,-4,-3,-10,-6,-9,-9,-7,-13,-3,-2,1,-16,-15,0,-5,-5,-3,-17,
        -7,7,-5,-11,-8,-2,-7,-2,0,-6,-6,2,-3,-12,-4,-2,-2,-5,-9,-10,-7,-3,-3,-17,
        -7,-5,8,2,-9,-6,-6,-7,0,-6,-12,0,-10,-12,-9,1,0,-17,-3,-10,6,-6,-3,-17,
        -4,-11,2,8,-14,-10,0,-2,-3,-11,-15,-7,-13,-15,-13,-5,-6,-16,-6,-5,7,0,-3,-17,
        -10,-8,-9,-14,11,-16,-15,-5,-7,-11,-9,-13,-14,0,-12,-1,-6,-2,0,-8,-10,-16,-5,-17,
        -11,-2,-6,-10,-16,8,-2,-10,0,-12,-4,0,-8,-12,-1,-9,-8,-14,-9,-13,-7,6,-4,-17,
        -4,-7,-6,0,-15,-2,7,-1,-9,-12,-15,-1,-10,-17,-13,-11,-8,-15,-12,-5,0,6,-4,-17,
        -3,-2,-7,-2,-5,-10,-1,7,-10,-11,-14,-6,-12,-9,-11,-1,-7,-5,-14,-5,-4,-3,-4,-17,
        -10,0,0,-3,-7,0,-9,-10,10,-10,-4,-5,-10,-6,-3,-6,-6,-11,2,-14,-1,-2,-3,-17,
        -6,-6,-6,-11,-11,-12,-12,-11,-10,7,0,-7,0,-2,-10,-4,0,-14,-9,2,-7,-12,-2,-17,
        -9,-6,-12,-15,-9,-4,-15,-14,-4,0,6,-10,0,0,-3,-5,-8,-6,-8,-4,-13,-6,-4,-17,
        -9,2,0,-7,-13,0,-1,-6,-5,-7,-10,7,-4,-14,-9,-5,-1,-12,-13,-9,-1,-1,-2,-17,
        -7,-3,-10,-13,-14,-8,-10,-12,-10,0,0,-4,10,-7,-11,-9,-1,-11,-15,0,-11,-9,-3,-17,
        -13,-12,-12,-15,0,-12,-17,-9,-6,-2,0,-14,-7,10,-11,-5,-10,-5,1,-5,-13,-14,-3,-17,
        -3,-4,-9,-13,-12,-1,-13,-11,-3,-10,-3,-9,-11,-11,8,-1,-3,-13,-11,-12,-10,-3,-5,-17,
        -2,-2,1,-5,-1,-9,-11,-1,-6,-4,-5,-5,-9,-5,-1,8,0,-12,-6,-9,0,-10,-3,-17,
        1,-2,0,-6,-6,-8,-8,-7,-6,0,-8,-1,-1,-10,-3,0,7,-16,-10,-4,-2,-8,-2,-17,
        -16,-5,-17,-16,-2,-14,-15,-5,-11,-14,-6,-12,-11,-5,-13,-12,-16,10,-4,-16,-16,-14,-8,-17,
        -15,-9,-3,-6,0,-9,-12,-14,2,-9,-8,-13,-15,1,-11,-6,-10,-4,10,-12,-4,-10,-4,-17,
        0,-10,-10,-5,-8,-13,-5,-5,-14,2,-4,-9,0,-5,-12,-9,-4,-16,-12,7,-7,-7,-3,-17,
        -5,-7,6,7,-10,-7,0,-4,-1,-7,-13,-1,-11,-13,-10,0,-2,-16,-4,-7,7,-2,-4,-17,
        -5,-3,-6,0,-16,6,6,-3,-2,-12,-6,-1,-9,-14,-3,-10,-8,-14,-10,-7,-2,6,-4,-17,
        -3,-3,-3,-3,-5,-4,-4,-4,-3,-2,-4,-2,-3,-3,-5,-3,-2,-8,-4,-3,-4,-4,-3,-17,
        -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,1
    ]).unwrap();
}

#[inline]
fn lookup(a: u8) -> usize {
    if a == b'Y' {
        23 as usize
    } else if a == b'Z' {
        24 as usize
    } else if a == b'X' {
        25 as usize
    } else if a == b'*' {
        26 as usize
    } else {
        (a - 65) as usize
    }
}

pub fn hiv25(a: u8, b: u8) -> i32 {
    let a = lookup(a);
    let b = lookup(b);

    MAT[(a, b)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hiv25() {
        let score1 = hiv25(b'A', b'A');
        assert_eq!(score1, 7);
    }
}
