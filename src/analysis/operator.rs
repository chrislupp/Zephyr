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

//! Operators evaluated at stored targets, with reflected source images on demand.
use super::Work;
use crate::core::{MatX, Scalar, Vec3};
use crate::solver::{VortexLattice, horseshoe_vortex_velocity};

pub(crate) fn reflect(v: Vec3) -> Vec3 {
    Vec3::new(v.x, -v.y, v.z)
}

pub(crate) fn core(mesh: &VortexLattice, row: usize, column: usize) -> Scalar {
    let source = &mesh.vortices()[column];
    let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
    if mesh.vortices()[row].component_index == source.component_index {
        0.0001 * width
    } else {
        (0.0001 * source.local_chord).max(0.0001 * width)
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble(
    mesh: &VortexLattice,
    mach: Scalar,
    normals: &[Vec3],
    half: bool,
    aic: bool,
    velocity: bool,
    work: &mut Work,
) -> (MatX, MatX) {
    let n = mesh.vortex_count();
    let pg = (1.0 - mach * mach).sqrt();
    let mut matrix = if aic {
        MatX::zeros(n, n)
    } else {
        MatX::zeros(0, 0)
    };
    let mut velocities = if velocity {
        MatX::zeros(3 * n, n)
    } else {
        MatX::zeros(0, 0)
    };
    for (j, source) in mesh.vortices().iter().enumerate() {
        for (i, target) in mesh.vortices().iter().enumerate() {
            let radius = core(mesh, i, j);
            let mut v = horseshoe_vortex_velocity(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                radius,
                i != j,
            );
            if half {
                // This is an explicit duplicate's influence: its bound leg is
                // included even for a coincident image midpoint, as in full mode.
                v += horseshoe_vortex_velocity(
                    &target.control_point,
                    &reflect(source.bound_vortex_right),
                    &reflect(source.bound_vortex_left),
                    pg,
                    radius,
                    true,
                );
            }
            if aic {
                matrix[(i, j)] = v.x * normals[i].x + v.y * normals[i].y + v.z * normals[i].z;
            }
            if velocity {
                for k in 0..3 {
                    velocities[(3 * i + k, j)] = v[k];
                }
            }
        }
    }
    work.influence_evaluations += n * n * if half { 2 } else { 1 };
    (matrix, velocities)
}
