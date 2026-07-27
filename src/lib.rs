/*
    Project Pegasus

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

#![forbid(unsafe_code)]
#![doc = include_str!("../README.md")]

//! # Conventions
//!
//! Pegasus uses aircraft body axes: X aft, Y starboard, and Z up. Flight
//! angles and rates are radians. Section incidence and surface angle offsets
//! retain AVL's historical degree convention.

pub mod aero;
pub mod core;
pub mod dynamics;
pub mod error;
pub mod geometry;
pub mod gradients;
pub mod io;
pub mod solver;
pub mod trim;
pub mod utils;
pub mod visualization;

pub use error::{ParseLocation, PegasusError, Result};
