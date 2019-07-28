// Copyright 2019, Todd Stellanova
// License: see LICENSE file


/// Rust implementation of the Arc* algorithm 1 described in:
/// "Asynchronous Corner Detection and Tracking for Event Cameras in Real Time", Alzugaray & M. Chli,
/// IEEE Robotics and Automation letters 2018  accessed via:
///[ETH Zurich archive](https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/277131/RAL2018-camera-ready.pdf?sequence=1&isAllowed=y)
///
/// Original Arc* algorithm described in pseudocode roughly as:
/// ```ignore
///    let mut Anew  = NewestElement(c3_vals);
///    let mut ECW = NextElementCW(Anew, c3_vals);
///    let mut ECCW = NextElementCCW(Anew, c3_vals);
///
///   while ECW != ECCW {
///        if ECW > ECCW {
///            if OldestElement(Anew) <= ECW || Length(Anew) < Lmin {
///                ExpandUntilElement(Anew, ECW);
///            }
///            ECW = NextElementCW(ECW, c3_vals);
///        } else {
///            if OldestElement(Anew) <= ECCW || Length(Anew) < Lmin {
///                ExpandUntilElement(Anew, ECCW);
///            }
///            ECCW = NextElementCCW(ECCW, c3_vals);
///        }
///    }
///    let newest_segment_size = Length(Anew);
///    if (Lmin <= newest_segment_size && newest_segment_size <= Lmax) ||  Lmin <= Length(C \ Anew ) <= Lmax {
///        return  true;
///    }
/// ```
/// This algorithm works with the Surface of Active Events (SAE) which is a matrix of
/// timestamps (one per pixel), indicating when a change event (rising or falling above or
/// below the detection threshold) most recently triggered at a particular pixel.

use arrayvec::ArrayVec;
use crate::sae_types::*;


const CIRCLE3_DIM: usize = 16;
/// pixel offsets of radius 3 circle surrounding point of interest
const CIRCLE3_GEN: [[i32; 2] ; CIRCLE3_DIM] = [
    [0, 3], [1, 3], [2, 2], [3, 1],
    [3, 0], [3, -1], [2, -2], [1, -3],
    [0, -3], [-1, -3], [-2, -2], [-3, -1],
    [-3, 0], [-3, 1], [-2, 2], [-1, 3]
];
const CIRCLE3_MIN_ARC_LEN:usize = 3;
const CIRCLE3_MAX_ARC_LEN:usize = 6;
type Circle3Vals = ArrayVec<[SaeTime;CIRCLE3_DIM]>;

const CIRCLE4_DIM: usize = 20;
/// pixel offsets of radius 4 circle surrounding point of interest
const CIRCLE4_GEN: [[i32; 2] ; CIRCLE4_DIM]  = [
    [0, 4], [1, 4], [2, 3], [3, 2],
    [4, 1], [4, 0], [4, -1], [3, -2],
    [2, -3], [1, -4], [0, -4], [-1, -4],
    [-2, -3], [-3, -2], [-4, -1], [-4, 0],
    [-4, 1], [-3, 2], [-2, 3], [-1, 4]
];
const CIRCLE4_MIN_ARC_LEN:usize = 4;
const CIRCLE4_MAX_ARC_LEN:usize = 8;
type Circle4Vals = ArrayVec<[SaeTime;CIRCLE4_DIM]>;


/// Number of pixels inset from all borders where we can start evaluating corners
const BORDER_INSET: usize = 4;

/// Get array of SAE values from the C3 circle surrounding the given point
fn c3_vals_for_point(sae_pol: &SaeMatrix, row: usize, col: usize) -> Circle3Vals {
    let mut res = Circle3Vals::new();

    let irow = row as i32;
    let icol = col as i32;

    for item in CIRCLE3_GEN.iter() {
        let a = (item[0] + irow) as usize;
        let b = (item[1] + icol) as usize;
        res.push(sae_pol[(a, b)] );
    }

    res
}

/// Get array of SAE values from the C4circle surrounding the given point
fn c4_vals_for_point(sae_pol: &SaeMatrix, row: usize, col: usize) -> Circle4Vals {
    let mut res = Circle4Vals::new();

    let irow = row as i32;
    let icol = col as i32;

    for item in CIRCLE4_GEN.iter() {
        let a = (item[0] + irow) as usize;
        let b = (item[1] + icol) as usize;
        res.push(sae_pol[(a, b)] );
    }

    res
}



/// Find the freshest timestamp in the given circle
fn find_freshest_in_circle(circle_vals: &[SaeTime]) -> (usize, SaeTime) {
    let mut newest_idx = 0;
    let mut newest_val: SaeTime = 0;
    //find the newest val in the circle
    for i in 0..circle_vals.len() {
        let val = circle_vals[i];
        if val > newest_val {
            newest_val = val;
            newest_idx = i;
        }
    }

    (newest_idx, newest_val)
}

/// returns the size of the arc segment containing the freshest SAE timestamps
fn arcstar_expand(circle_vals: &[SaeTime], circle_dim: usize, min_arc_size: usize,  newest_idx: usize)  -> usize {

    let mut cw_idx:usize = (newest_idx + 1) % circle_dim;
    let mut ccw_idx:usize = (newest_idx + (circle_dim-1)) % circle_dim;

    let mut arc_cw_val = circle_vals[cw_idx];
    let mut arc_ccw_val = circle_vals[ccw_idx];
    let mut arc_cw_oldest = arc_cw_val;
    let mut arc_ccw_oldest = arc_ccw_val;
    let mut segment_oldest =  SaeTime::max_value();

    //Expand beginning with pixels immediately neighboring newest_idx
    for _iteration in 1..min_arc_size {
        // Pick CW/CCW expansion based on which next circle item has freshest timestamp
        if arc_cw_val > arc_ccw_val {
            // CW arc has freshest value: include arc in new segment
            if arc_cw_oldest < segment_oldest {
                segment_oldest = arc_cw_oldest;
            }
            // Expand arc cw
            cw_idx = ( cw_idx + 1 ) % circle_dim;
            arc_cw_val = circle_vals[cw_idx];
            if arc_cw_val < arc_cw_oldest {
                // Update oldest item in the arc
                arc_cw_oldest = arc_cw_val;
            }
        }
        else {
            // CCW arc has freshest value: include arc in new segment
            if arc_ccw_oldest < segment_oldest {
                segment_oldest = arc_ccw_oldest;
            }
            // Expand arc ccw
            ccw_idx = (ccw_idx + (circle_dim - 1)) % circle_dim;
            arc_ccw_val = circle_vals[ccw_idx];
            if arc_ccw_val < arc_ccw_oldest {
                // Update oldest item in the arc
                arc_ccw_oldest = arc_ccw_val;
            }
        }
    }

    // this is the arc length of the arc containing the freshest elements in the circle
    //TODO check this assumption
    let mut freshest_arc_size: usize = min_arc_size;

    // Continue expansion, looking at freshest values
    for iteration in min_arc_size..circle_dim {
        // Pick CW/CCW expansion based on which next circle item has freshest timestamp
        if arc_cw_val > arc_ccw_val {
            // CW arc has the freshest value: include arc in freshest segment
            if arc_cw_val >=  segment_oldest {
                freshest_arc_size = iteration + 1;
                if arc_cw_oldest < segment_oldest {
                    segment_oldest = arc_cw_oldest;
                }
            }
            // Expand arc clockwise
            cw_idx = ( cw_idx + 1) % circle_dim;
            arc_cw_val = circle_vals[cw_idx];
            if arc_cw_val < arc_cw_oldest {
                // Update oldest item in the arc
                arc_cw_oldest = arc_cw_val;
            }
        }
        else {
            // CCW arc has the freshest value: include arc in freshest segment
            if arc_ccw_val >=  segment_oldest {
                freshest_arc_size = iteration + 1;
                if arc_ccw_oldest < segment_oldest {
                    segment_oldest = arc_ccw_oldest;
                }
            }
            // Expand arc counter-clockwise
            ccw_idx = (ccw_idx + (circle_dim - 1) ) % circle_dim;
            arc_ccw_val = circle_vals[ccw_idx];
            if arc_ccw_val < arc_ccw_oldest {
                // Update oldest item in the arc
                arc_ccw_oldest = arc_ccw_val;
            }
        }
    }

    freshest_arc_size
}

/// returns whether the given point in updated SAE is a corner
fn arcstar_check_for_point(sae_pol: &SaeMatrix, evt: &mut SaeEvent) -> bool {
    let row = evt.row as usize;
    let col = evt.col as usize;

    let c3_vals:Circle3Vals = c3_vals_for_point(sae_pol, row, col);
    let c3_vals_slice = c3_vals.as_slice();
    let (freshest_c3_idx, freshest_c3_val) = find_freshest_in_circle(c3_vals_slice);
    let freshest_c3_segment_size = arcstar_expand(c3_vals_slice, CIRCLE3_DIM, CIRCLE3_MIN_ARC_LEN, freshest_c3_idx);

    let mut arc_valid =
        (freshest_c3_segment_size <= CIRCLE3_MAX_ARC_LEN) ||
            ((freshest_c3_segment_size >= (CIRCLE3_DIM - CIRCLE3_MAX_ARC_LEN)) &&
                (freshest_c3_segment_size <= (CIRCLE3_DIM - CIRCLE3_MIN_ARC_LEN) ));

    if arc_valid {
        let c4_vals:Circle4Vals = c4_vals_for_point(sae_pol, row, col);
        let c4_vals_slice = c4_vals.as_slice();

        let (freshest_c4_idx, freshest_c4_val) = find_freshest_in_circle(c4_vals_slice);
        let freshest_c4_segment_size = arcstar_expand(c4_vals_slice, CIRCLE4_DIM, CIRCLE4_MIN_ARC_LEN, freshest_c4_idx);
        arc_valid =
            (freshest_c4_segment_size <= CIRCLE4_MAX_ARC_LEN) ||
                ((freshest_c4_segment_size >= (CIRCLE4_DIM - CIRCLE4_MAX_ARC_LEN)) &&
                    (freshest_c4_segment_size <= (CIRCLE4_DIM - CIRCLE4_MIN_ARC_LEN) ));

        if arc_valid {
            //this is where we calculate the descriptor "fingerprint" for an event,
            //based on the shape of the surrounding SAE
            let freshest_seg_val:f32 = (freshest_c3_val.max(freshest_c4_val)) as f32;
            let mut desc_idx = 0;
            let mut norm_descriptor:NormDescriptor = [0.0; NORM_DESCRIPTOR_LEN];
            //iterate around C3 starting from maximum index
            for c3_idx in 0..c3_vals_slice.len() {
                let true_idx = (c3_idx + freshest_c3_idx) % c3_vals_slice.len();
                let val = c3_vals_slice[true_idx];
                let norm: f32 = 1.0f32 - (freshest_seg_val - (val as f32))/freshest_seg_val;
                norm_descriptor[desc_idx] = norm;
                desc_idx +=1;
            }
            //iterate around C4 starting from maximum index
            for c4_idx in 0..c4_vals_slice.len() {
                let true_idx = (c4_idx + freshest_c4_idx) % c4_vals_slice.len();
                let val = c4_vals_slice[true_idx];
                let norm: f32 = 1.0f32 - (freshest_seg_val - (val as f32))/freshest_seg_val;
                norm_descriptor[desc_idx] = norm;
                desc_idx +=1;
            }

            evt.norm_descriptor = Some(Box::new(norm_descriptor));
        }
    }

    arc_valid
}

fn arcstar_is_event_corner(sae_pol: &SaeMatrix, evt: &mut SaeEvent) -> bool {
    let row = evt.row as usize;
    let col = evt.col as usize;

    //filter out events too close to SAE border
    let (nrows, ncols) = sae_pol.shape();
    if (col < BORDER_INSET) || (col >= (ncols - BORDER_INSET)) ||
        (row < BORDER_INSET) || (row >= (nrows - BORDER_INSET))  {
        //println!("shape: ({}, {}) border: (row {} , col {})" , nrows, ncols, row, col);
        return false;
    }

    arcstar_check_for_point(sae_pol, evt)
}


/// Detect whether the input event is a corner, and compute descriptor if so:
/// returns a modified event with computed descriptor, if it's a corner.
pub fn detect_and_compute_one(sae_pol: &SaeMatrix, evt: &SaeEvent) -> Option<SaeEvent> {
    let mut out_evt: SaeEvent = evt.clone();

    match arcstar_is_event_corner(sae_pol, &mut out_evt) {
        true => Some(out_evt),
        false => None
    }
}



#[cfg(test)]
mod tests {

    use super::*;

    type StaticSaeArray = [[SaeTime; 9] ; 9];

    const SAE_ALL_RAYS: StaticSaeArray = [
        [0, 0, 0, 0, 7, 0, 0, 0, 0 ],
        [0, 7, 0, 0, 7, 0, 0, 7, 0 ],
        [0, 0, 7, 0, 0, 0, 7, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 0, 0, 9, 0, 0, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 7, 0, 0, 0, 7, 0, 0 ],
        [0, 7, 0, 0, 7, 0, 0, 7, 0 ],
        [0, 0, 0, 0, 7, 0, 0, 0, 0 ]
    ];

    const SAE_BLANK: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_OUTSIDE_CORNER_NE: StaticSaeArray = [
        [0, 0, 0, 0, 80, 79, 78, 77, 76 ],
        [0, 0, 0, 0, 85, 84, 83, 82, 81 ],
        [0, 0, 0, 0, 90, 89, 88, 87, 86 ],
        [0, 0, 0, 0, 95, 94, 93, 92, 91 ],
        [0, 0, 0, 0, 100, 99, 98, 97, 96 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_OUTSIDE_CORNER_SE: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ]
    ];

    const SAE_OUTSIDE_CORNER_SW: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ]
    ];


    const SAE_OUTSIDE_CORNER_NW: StaticSaeArray = [
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];


    const SAE_OUTSIDE_CORNER_SSE: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 7, 0 ],
        [0, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 7, 7, 7, 7, 7, 7, 7 ]
    ];
    const SAE_INSIDE_CORNER_NE: StaticSaeArray = [
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ]
    ];

    const SAE_INSIDE_CORNER_NW: StaticSaeArray = [
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 9, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ]
    ];

    const SAE_INSIDE_CORNER_SE: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ]
    ];

    const SAE_INSIDE_CORNER_SW: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 9, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ]
    ];

    const SAE_INSIDE_CORNER_N: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 0, 0, 0, 0, 0, 0, 0, 7 ],
        [7, 7, 0, 0, 0, 0, 0, 7, 7 ],
        [7, 7, 7, 0, 0, 0, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ]
    ];

    const SAE_INSIDE_CORNER_S: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 7, 7, 7, 7 ],
        [7, 7, 7, 0, 0, 0, 7, 7, 7 ],
        [7, 7, 0, 0, 0, 0, 0, 7, 7 ],
        [7, 0, 0, 0, 0, 0, 0, 0, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_INSIDE_CORNER_E: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 0 ]
    ];

    const SAE_INSIDE_CORNER_W: StaticSaeArray = [
        [0, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 9, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 7, 7, 7, 7, 7, 7, 7, 7 ]
    ];

    const SAE_OUTSIDE_CORNER_N: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 7, 7, 7, 7, 7, 7, 7, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    pub const SAE_OUTSIDE_CORNER_S: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 7, 7, 7, 7, 7, 7, 7, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ]
    ];

    const SAE_OUTSIDE_CORNER_E: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 9, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 7 ]
    ];

    const SAE_OUTSIDE_CORNER_W: StaticSaeArray = [
        [7, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_BAR_VERT_THICK: StaticSaeArray = [
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 0, 0, 0, 0 ]
    ];

    const SAE_BAR_VERT_THIN: StaticSaeArray = [
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 9, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 0, 0, 0, 0, 0 ]
    ];

    const SAE_CENTER_BAR_VERT_THICK: StaticSaeArray = [
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 9, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
        [0, 0, 7, 7, 7, 7, 7, 0, 0 ],
    ];


    const SAE_CENTER_BAR_VERT_THIN: StaticSaeArray = [
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 9, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
        [0, 0, 0, 7, 7, 7, 0, 0, 0 ],
    ];

    const SAE_DIAG_BAR_VERT_THIN: StaticSaeArray = [
        [0, 0, 0, 0, 0, 7, 7, 7, 0 ],
        [0, 0, 0, 0, 0, 7, 7, 7, 0 ],
        [0, 0, 0, 0, 7, 7, 7, 0, 0 ],
        [0, 0, 0, 0, 7, 7, 7, 0, 0 ],
        [0, 0, 0, 7, 9, 7, 0, 0, 0 ],
        [0, 0, 7, 7, 7, 0, 0, 0, 0 ],
        [0, 0, 7, 7, 7, 0, 0, 0, 0 ],
        [0, 7, 7, 7, 0, 0, 0, 0, 0 ],
        [0, 7, 7, 7, 0, 0, 0, 0, 0 ],
    ];

    const SAE_DIAG_BAR_NE_THIN: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 1 ],
        [0, 0, 0, 0, 0, 0, 0, 3, 0 ],
        [0, 0, 0, 0, 0, 0, 5, 0, 0 ],
        [0, 0, 0, 0, 0, 7, 0, 0, 0 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 8, 0, 0, 0, 0, 0 ],
        [0, 0, 6, 0, 0, 0, 0, 0, 0 ],
        [0, 4, 0, 0, 0, 0, 0, 0, 0 ],
        [2, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_BAR_HORIZ_THIN: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 9, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_BAR_HORIZ_THICK: StaticSaeArray = [
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_CENTER_BAR_HORIZ_THIN: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    const SAE_CENTER_BAR_HORIZ_THICK: StaticSaeArray = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 9, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [7, 7, 7, 7, 7, 7, 7, 7, 7 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ],
        [0, 0, 0, 0, 0, 0, 0, 0, 0 ]
    ];

    fn generate_test_event() -> SaeEvent {
        SaeEvent {
            row: 4,
            col: 4,
            polarity: 0,
            timestamp: 0,
            norm_descriptor: Some(Box::new([666.0f32; NORM_DESCRIPTOR_LEN])),
        }
    }

    fn init_matrix_from_static_sae_array(input: &StaticSaeArray) -> SaeMatrix {
        let mut sae_pol = SaeMatrix::zeros( 9, 9);

        for row in 0..9 {
            for col in 0..9 {
                sae_pol[(row, col)] = input[row][col];
            }
        }
        println!("{}", sae_pol);
        sae_pol
    }


    #[test]
    fn test_is_event_corner_blank() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_BLANK);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_outside_corner_ne() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_NE);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_outside_corner_se() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_SE);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_outside_corner_sw() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_SW);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_outside_corner_nw() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_NW);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_outside_corner_sse() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_SSE);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_n() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_N);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_s() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_S);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_e() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_E);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_w() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_W);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_ne() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_NE);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_inside_corner_nw() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_NW);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_inside_corner_se() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_SE);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_inside_corner_sw() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_INSIDE_CORNER_SW);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_outside_corner_north() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_N);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_outside_corner_south() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_S);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_outside_corner_east() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_E);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_outside_corner_west() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_OUTSIDE_CORNER_W);
        let mut evt = generate_test_event();
        assert_eq!(true, arcstar_is_event_corner(&sae_pol, &mut evt));    }

    #[test]
    fn test_is_event_bar_vertical() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_BAR_VERT_THICK);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_BAR_VERT_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_CENTER_BAR_VERT_THICK);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_CENTER_BAR_VERT_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));
    }


    #[test]
    fn test_is_event_diag_bar() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_DIAG_BAR_VERT_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_DIAG_BAR_NE_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_bar_horizontal() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_BAR_HORIZ_THICK);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_BAR_HORIZ_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_CENTER_BAR_HORIZ_THICK);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));

        let sae_pol = init_matrix_from_static_sae_array(&SAE_CENTER_BAR_HORIZ_THIN);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));
    }

    #[test]
    fn test_is_event_corner_all_rays() {
        let sae_pol = init_matrix_from_static_sae_array(&SAE_ALL_RAYS);
        let mut evt = generate_test_event();
        assert_eq!(false, arcstar_is_event_corner(&sae_pol, &mut evt));
    }


}
