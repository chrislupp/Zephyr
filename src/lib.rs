/*
    Project Zephyr

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-2198.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/

#![forbid(unsafe_code)]
#![doc = include_str!("../README.md")]

//! # Conventions
//!
//! Zephyr uses aircraft body axes: X aft, Y starboard, and Z up. Flight
//! angles and rates are radians. Section incidence and surface angle offsets
//! retain AVL's historical degree convention.

pub mod aero;
pub mod analysis;
pub mod core;
pub mod dynamics;
pub mod error;
pub mod geometry;
pub mod gradients;
pub mod io;
#[cfg(feature = "philote")]
pub mod philote;
pub mod solver;
pub mod trim;
pub mod utils;
pub mod visualization;

pub use error::{ParseLocation, Result, ZephyrError};
