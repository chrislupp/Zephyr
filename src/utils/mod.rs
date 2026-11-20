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

//! Mathematical helpers and interpolation.

mod cubic_spline;
mod math;

pub use cubic_spline::{CubicSpline, SplineBoundaryCondition};
pub use math::{
    blended_spacing, body_to_stability_rotation, cosine_spacing, prandtl_glauert_factor,
    rotate_about_x, rotate_about_y, rotate_about_z, sine_spacing, spacing,
    trailing_vortex_velocity, uniform_spacing, vortex_segment_velocity, wind_to_body_rotation,
};
