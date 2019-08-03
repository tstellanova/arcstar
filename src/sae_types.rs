// Copyright 2019, Todd Stellanova
// License: see LICENSE file

use nalgebra::{DMatrix};
use std::fmt;

/// Defining types for Surface of Active Events (SAE), commonly used for Event Cameras (DVS &c.).
/// The SAE is a matrix of  timestamps (one per pixel), indicating when a change event
/// (rising or falling above or below the detection threshold)
/// most recently triggered at a particular  pixel.

/// The type used to store timestamps in the SAE
pub type SaeTime = u32;
/// Type used to store a Surface of Active Events
pub type SaeMatrix = DMatrix<SaeTime>;


pub const NORM_DESCRIPTOR_LEN: usize = 36;
/// a crude feature descriptor allowing limited number of comparison points
pub type NormDescriptor = [f32; NORM_DESCRIPTOR_LEN];


/// The main change event struct
#[derive(Clone)]
pub struct SaeEvent {
  pub row: u16,
  pub col: u16,
  pub polarity: u8,
  pub timestamp: SaeTime,
  pub norm_descriptor: Option<Box<NormDescriptor>>,
}

impl fmt::Debug for SaeEvent {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    let mut avg_desc = 0.0;
    if self.norm_descriptor.is_some() {
      let values = self.norm_descriptor.clone().unwrap();
      let total:f32 = values.iter().sum();
      avg_desc = total / (NORM_DESCRIPTOR_LEN as f32);
    }
    write!(f, "SaeEvent {{ row: {}, col: {} time: {} pol: {} avg_desc: {} }}",
           self.row, self.col, self.timestamp, self.polarity, avg_desc)
  }
}

impl PartialEq for SaeEvent {
  fn eq(&self, other: &SaeEvent) -> bool {
    self.row == other.row &&
      self.col == other.col &&
      self.polarity == other.polarity &&
      self.timestamp == other.timestamp

    // TODO implement PartialEq for norm_descriptor ?
  }
}


impl std::default::Default for SaeEvent {
  fn default() -> Self {
    SaeEvent {
      row: 0,
      col: 0,
      polarity: 0,
      timestamp: 0,
      norm_descriptor: None,
    }
  }
}

impl SaeEvent {
  pub fn new() -> Self {
    Self::default()
  }

  /// square of the euclidean distance between two events
  pub fn spatial_dist_2(&self, other: &Self) -> u32 {
    let drow: u32 = (self.row.max(other.row) - self.row.min(other.row)) as u32;
    let dcol: u32 = (self.col.max(other.col) - self.col.min(other.col)) as u32;

    drow*drow + dcol*dcol
  }

  /// manhattan distance or rectilinear distance
  pub fn spatial_rl_dist(&self, other: &Self) -> u32 {
    let drow: u32 = (self.row.max(other.row) - self.row.min(other.row)) as u32;
    let dcol: u32 = (self.col.max(other.col) - self.col.min(other.col)) as u32;

    drow + dcol
  }

  pub fn likeness(&self, b: &SaeEvent) -> f32 {
    if self.norm_descriptor.is_none() || b.norm_descriptor.is_none() {
      return 0.0;
    }

    let b_desc = match b.norm_descriptor {
      Some(ref p) => p,
      _ => unreachable!()
    };

    let a_desc = match self.norm_descriptor {
      Some(ref p) => p,
      _ => unreachable!()
    };

    let mut da_total:f32 = 0.0;
    let mut db_total:f32 = 0.0;
    let mut min_total:f32 = 0.0;

    for i in 0..a_desc.len() {
      let da = a_desc[i];
      let db = b_desc[i];
      da_total += da;
      db_total += db;
      min_total += da.min(db);
    }

    let max_total = da_total.max(db_total);

    //likeness is 0..1
    let likeness: f32 = min_total/max_total;
    //println!("likeness: {}", likeness);

    likeness
  }
}



#[cfg(test)]
mod tests {
  use super::*;
  use assert_approx_eq::assert_approx_eq;

  #[test]
  fn test_spatial_dist() {
    let mut evt_a:SaeEvent = SaeEvent::new();
    let mut evt_b:SaeEvent = SaeEvent::new();
    let dist2 = evt_a.spatial_dist_2(&evt_b);
    let rl_dist = evt_a.spatial_rl_dist(&evt_b);
    assert_eq!(rl_dist, 0);
    assert_eq!(dist2, 0);

    evt_a.col = 3;
    evt_b.col = 5;

    let dist2 = evt_a.spatial_dist_2(&evt_b);
    let rl_dist = evt_a.spatial_rl_dist(&evt_b);
    assert_eq!(rl_dist, 2);
    assert_eq!(dist2, 4);

    evt_a.row = 3;
    evt_b.row = 5;

    let dist2 = evt_a.spatial_dist_2(&evt_b);
    let rl_dist = evt_a.spatial_rl_dist(&evt_b);
    assert_eq!(rl_dist, 4);
    assert_eq!(dist2, 8);

  }

  #[test]
  fn test_event_likeness() {
    let evt_a = SaeEvent {
      row: 0,
      col: 0,
      polarity: 0,
      timestamp: 0,
      norm_descriptor: Some(Box::new([1.0; NORM_DESCRIPTOR_LEN])),
    };

    let mut evt_b = SaeEvent {
      row: 0,
      col: 0,
      polarity: 0,
      timestamp: 0,
      norm_descriptor: Some(Box::new([1.0; NORM_DESCRIPTOR_LEN])),
    };

    let likeness = evt_a.likeness(&evt_b);
    assert_approx_eq!(likeness, 1.0);

    evt_b.norm_descriptor = Some(Box::new([0.0; NORM_DESCRIPTOR_LEN]));
    let likeness = evt_a.likeness(&evt_b);
    assert_approx_eq!(likeness, 0.0);

    evt_b.norm_descriptor = Some(Box::new([0.5; NORM_DESCRIPTOR_LEN]));
    let likeness = evt_a.likeness(&evt_b);
    assert_approx_eq!(likeness, 0.5);

    evt_b.norm_descriptor = None;
    let likeness = evt_a.likeness(&evt_b);
    assert_approx_eq!(likeness, 0.0);
  }

}