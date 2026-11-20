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

//! Analytical compressibility derivatives of the velocity influence matrix.

use crate::core::{EPSILON, MatX, SymmetryType, Vec3};
use crate::solver::VortexLatticeSolver;

use super::horseshoe_velocity_pg_derivative;

pub(crate) fn velocity_influence_mach_derivative(solver: &VortexLatticeSolver) -> MatX {
    let mesh = solver.mesh();
    let n = mesh.vortex_count();
    let mach = solver.current_mach();
    if !(0.0..1.0).contains(&mach) || mach <= EPSILON {
        return MatX::zeros(3 * n, n);
    }
    let pg = (1.0 - mach * mach).sqrt();
    let pg_direction = -mach / pg;
    let symmetry = solver.symmetry();
    let mut matrix = MatX::zeros(3 * n, n);
    for row in 0..n {
        let target = &mesh.vortices()[row];
        for column in 0..n {
            let source = &mesh.vortices()[column];
            let width = (source.bound_vortex_right - source.bound_vortex_left).norm();
            let same_component = target.component_index == source.component_index;
            let core = if same_component {
                0.0001 * width
            } else {
                (0.0001 * source.local_chord).max(0.0001 * width)
            };
            let mut derivative = horseshoe_velocity_pg_derivative(
                &target.control_point,
                &source.bound_vortex_left,
                &source.bound_vortex_right,
                pg,
                core,
                row != column,
            );
            if symmetry.y_symmetry != SymmetryType::None {
                let offset = 2.0 * symmetry.y_symmetry_plane;
                let left = Vec3::new(
                    source.bound_vortex_right.x,
                    offset - source.bound_vortex_right.y,
                    source.bound_vortex_right.z,
                );
                let right = Vec3::new(
                    source.bound_vortex_left.x,
                    offset - source.bound_vortex_left.y,
                    source.bound_vortex_left.z,
                );
                let include_bound = !(symmetry.y_symmetry == SymmetryType::Symmetric
                    && (target.control_point - 0.5 * (left + right)).norm_squared() < EPSILON);
                derivative += symmetry.y_symmetry.factor()
                    * horseshoe_velocity_pg_derivative(
                        &target.control_point,
                        &left,
                        &right,
                        pg,
                        core,
                        include_bound,
                    );
            }
            if symmetry.z_symmetry != SymmetryType::None {
                let offset = 2.0 * symmetry.z_symmetry_plane;
                let left = Vec3::new(
                    source.bound_vortex_right.x,
                    source.bound_vortex_right.y,
                    offset - source.bound_vortex_right.z,
                );
                let right = Vec3::new(
                    source.bound_vortex_left.x,
                    source.bound_vortex_left.y,
                    offset - source.bound_vortex_left.z,
                );
                derivative += symmetry.z_symmetry.factor()
                    * horseshoe_velocity_pg_derivative(
                        &target.control_point,
                        &left,
                        &right,
                        pg,
                        core,
                        true,
                    );
                if symmetry.y_symmetry != SymmetryType::None {
                    let y_offset = 2.0 * symmetry.y_symmetry_plane;
                    let z_offset = 2.0 * symmetry.z_symmetry_plane;
                    let left = Vec3::new(
                        source.bound_vortex_left.x,
                        y_offset - source.bound_vortex_left.y,
                        z_offset - source.bound_vortex_left.z,
                    );
                    let right = Vec3::new(
                        source.bound_vortex_right.x,
                        y_offset - source.bound_vortex_right.y,
                        z_offset - source.bound_vortex_right.z,
                    );
                    derivative += symmetry.y_symmetry.factor()
                        * symmetry.z_symmetry.factor()
                        * horseshoe_velocity_pg_derivative(
                            &target.control_point,
                            &left,
                            &right,
                            pg,
                            core,
                            true,
                        );
                }
            }
            derivative *= pg_direction;
            matrix[(3 * row, column)] = derivative.x;
            matrix[(3 * row + 1, column)] = derivative.y;
            matrix[(3 * row + 2, column)] = derivative.z;
        }
    }
    matrix
}
