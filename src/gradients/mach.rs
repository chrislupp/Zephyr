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
//!
//! The derivative matrix `dV/dM` is `3N x N` and dense, but every consumer
//! multiplies it by the circulation straight away. The product is therefore
//! evaluated row by row without materialising the matrix: each target row
//! accumulates its `N` source contributions in column order, which is exactly
//! the order a column-major matrix-vector product adds them in, and no row reads
//! anything another row writes. That makes the rows free to run on any thread
//! in any order with a bit-identical result.

use crate::core::{EPSILON, Scalar, SymmetryType, Vec3, VecX};
use crate::solver::{SymmetryConfig, VortexLattice, VortexLatticeSolver};

use super::horseshoe_velocity_pg_derivative;

#[cfg(test)]
use crate::core::MatX;

/// Panel count below which the product stays on one thread.
///
/// Matches the assembly driver's threshold: the product is `O(N^2)` while
/// handing out rows is `O(N)`. Affects speed only -- both branches evaluate
/// every row with the same code.
#[cfg(feature = "parallel")]
const PARALLEL_MIN_PANELS: usize = 128;

/// `dV/dM * circulation`, the Mach derivative of the induced velocity at every
/// control point (`at_bound = false`) or bound-vortex centre (`at_bound = true`)
/// for fixed circulation, as a flat `3N` vector of x/y/z triples.
///
/// Identically zero, and computed without touching the lattice, at Mach zero
/// or outside the subsonic range the Prandtl-Glauert factor is defined on.
pub(crate) fn velocity_influence_mach_derivative_times(
    solver: &VortexLatticeSolver,
    at_bound: bool,
    circulation: &VecX,
) -> VecX {
    let n = solver.mesh().vortex_count();
    let Some(derivative) = MachDerivative::new(solver, at_bound) else {
        return VecX::zeros(3 * n);
    };
    let rows = drive_rows(n, |row| derivative.row_times(row, circulation));
    VecX::from_fn(3 * n, |index, _| rows[index / 3][index % 3])
}

impl VortexLatticeSolver {
    /// [`velocity_influence_mach_derivative_times`] for crate code outside
    /// `gradients`, whose `mach` module is private.
    pub(crate) fn velocity_mach_derivative_times(
        &self,
        at_bound: bool,
        circulation: &VecX,
    ) -> VecX {
        velocity_influence_mach_derivative_times(self, at_bound, circulation)
    }
}

/// The dense `3N x N` matrix whose product the function above evaluates.
#[cfg(test)]
pub(crate) fn velocity_influence_mach_derivative(
    solver: &VortexLatticeSolver,
    at_bound: bool,
) -> MatX {
    let n = solver.mesh().vortex_count();
    let mut matrix = MatX::zeros(3 * n, n);
    if let Some(derivative) = MachDerivative::new(solver, at_bound) {
        for row in 0..n {
            let point = derivative.field_point(row);
            for column in 0..n {
                let entry = derivative.entry(row, &point, column);
                matrix[(3 * row, column)] = entry.x;
                matrix[(3 * row + 1, column)] = entry.y;
                matrix[(3 * row + 2, column)] = entry.z;
            }
        }
    }
    matrix
}

/// Evaluate `row_value` for every row, in row order.
#[cfg(not(feature = "parallel"))]
fn drive_rows(count: usize, row_value: impl Fn(usize) -> Vec3 + Sync) -> Vec<Vec3> {
    (0..count).map(row_value).collect()
}

/// Evaluate `row_value` for every row, one rayon task per row. Each row is a
/// self-contained sequential sum and `collect` keeps row order, so the result is
/// the serial one bit for bit, whatever the thread count.
#[cfg(feature = "parallel")]
fn drive_rows(count: usize, row_value: impl Fn(usize) -> Vec3 + Sync) -> Vec<Vec3> {
    use rayon::prelude::*;

    if count < PARALLEL_MIN_PANELS {
        return (0..count).map(row_value).collect();
    }
    (0..count).into_par_iter().map(&row_value).collect()
}

/// Everything one entry of `dV/dM` reads, with the per-source core radii
/// hoisted out of the `N^2` loop.
struct MachDerivative<'a> {
    mesh: &'a VortexLattice,
    symmetry: SymmetryConfig,
    at_bound: bool,
    prandtl_glauert_factor: Scalar,
    /// `d(pg)/dM`, which turns the kernel's `d/d(pg)` into `d/dM`.
    pg_direction: Scalar,
    /// Core radius of each source panel seen from its own component.
    core_same_component: Vec<Scalar>,
    /// Core radius of each source panel seen from any other component.
    core_other_component: Vec<Scalar>,
}

impl<'a> MachDerivative<'a> {
    fn new(solver: &'a VortexLatticeSolver, at_bound: bool) -> Option<Self> {
        let mach = solver.current_mach();
        if !(0.0..1.0).contains(&mach) || mach <= EPSILON {
            return None;
        }
        let mesh = solver.mesh();
        let pg = (1.0 - mach * mach).sqrt();
        let (core_same_component, core_other_component) = mesh
            .vortices()
            .iter()
            .map(|source| {
                let width = (source.bound_vortex_right.y - source.bound_vortex_left.y)
                    .hypot(source.bound_vortex_right.z - source.bound_vortex_left.z);
                (
                    0.0001 * width,
                    (0.0001 * source.local_chord).max(0.0001 * width),
                )
            })
            .unzip();
        Some(Self {
            mesh,
            symmetry: solver.symmetry(),
            at_bound,
            prandtl_glauert_factor: pg,
            pg_direction: -mach / pg,
            core_same_component,
            core_other_component,
        })
    }

    fn field_point(&self, row: usize) -> Vec3 {
        let target = &self.mesh.vortices()[row];
        if self.at_bound {
            target.bound_vortex_center
        } else {
            target.control_point
        }
    }

    /// Row `row` of `dV/dM` times the circulation, summed in column order.
    fn row_times(&self, row: usize, circulation: &VecX) -> Vec3 {
        let point = self.field_point(row);
        let mut sum = Vec3::zeros();
        for (column, gamma) in circulation.iter().enumerate() {
            sum += *gamma * self.entry(row, &point, column);
        }
        sum
    }

    /// The `3 x 1` block of `dV/dM` coupling source `column` to target `row`,
    /// including the symmetry images the solver assembles.
    fn entry(&self, row: usize, point: &Vec3, column: usize) -> Vec3 {
        let pg = self.prandtl_glauert_factor;
        let symmetry = &self.symmetry;
        let target = &self.mesh.vortices()[row];
        let source = &self.mesh.vortices()[column];
        let core = if target.component_index == source.component_index {
            self.core_same_component[column]
        } else {
            self.core_other_component[column]
        };
        let pg_derivative = |left: &Vec3, right: &Vec3, include_bound: bool| {
            horseshoe_velocity_pg_derivative(point, left, right, pg, core, include_bound)
        };
        let mut derivative = pg_derivative(
            &source.bound_vortex_left,
            &source.bound_vortex_right,
            !self.at_bound || row != column,
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
                && (point - 0.5 * (left + right)).norm_squared() < EPSILON);
            derivative +=
                symmetry.y_symmetry.factor() * pg_derivative(&left, &right, include_bound);
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
            derivative += symmetry.z_symmetry.factor() * pg_derivative(&left, &right, true);
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
                    * pg_derivative(&left, &right, true);
            }
        }
        derivative * self.pg_direction
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Vec3;
    use crate::geometry::{Aircraft, Section, Surface};
    use crate::solver::FlightCondition;

    /// A wing and a tail on separate components, so both core-radius branches
    /// are taken.
    fn aircraft(spanwise: usize) -> Aircraft {
        let mut wing = Surface::new("Wing");
        wing.set_chordwise_panel_count(3);
        wing.set_spanwise_panel_count(spanwise);
        wing.set_component_index(0);
        wing.add_section(Section::new(Vec3::new(0.0, 0.1, 0.0), 1.0, 2.0));
        wing.add_section(Section::new(Vec3::new(0.3, 2.0, 0.2), 0.6, 0.0));
        let mut tail = Surface::new("Tail");
        tail.set_chordwise_panel_count(2);
        tail.set_spanwise_panel_count(3);
        tail.set_component_index(1);
        tail.add_section(Section::new(Vec3::new(3.0, 0.1, 0.4), 0.5, 0.0));
        tail.add_section(Section::new(Vec3::new(3.2, 0.9, 0.5), 0.3, 0.0));
        let mut aircraft = Aircraft::new("Mach derivative fixture");
        aircraft.add_surface(wing);
        aircraft.add_surface(tail);
        aircraft.set_reference_area(4.0);
        aircraft.set_reference_chord(1.0);
        aircraft.set_reference_span(4.0);
        aircraft
    }

    fn solver(mach: f64, symmetry: SymmetryConfig) -> VortexLatticeSolver {
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft(4)).unwrap();
        solver.set_symmetry(symmetry);
        solver
            .solve(&FlightCondition {
                alpha: 0.05,
                mach,
                ..FlightCondition::default()
            })
            .unwrap();
        solver
    }

    fn symmetries() -> [SymmetryConfig; 4] {
        [
            SymmetryConfig::default(),
            SymmetryConfig {
                y_symmetry: SymmetryType::Symmetric,
                ..SymmetryConfig::default()
            },
            SymmetryConfig {
                z_symmetry: SymmetryType::Antisymmetric,
                z_symmetry_plane: -0.7,
                ..SymmetryConfig::default()
            },
            SymmetryConfig {
                y_symmetry: SymmetryType::Antisymmetric,
                y_symmetry_plane: -0.05,
                z_symmetry: SymmetryType::Symmetric,
                z_symmetry_plane: -0.9,
            },
        ]
    }

    fn circulation(n: usize) -> VecX {
        VecX::from_fn(n, |i, _| 0.2 + 0.05 * (i as f64).sin())
    }

    #[test]
    fn the_product_is_the_dense_matrix_times_the_circulation() {
        for symmetry in symmetries() {
            let solver = solver(0.45, symmetry);
            let gamma = circulation(solver.mesh().vortex_count());
            for at_bound in [false, true] {
                let product = velocity_influence_mach_derivative_times(&solver, at_bound, &gamma);
                let dense = velocity_influence_mach_derivative(&solver, at_bound) * &gamma;
                let scale = dense.amax();
                assert!(scale > 0.0);
                assert!((product - dense).amax() <= 1.0e-14 * scale);
            }
        }
    }

    #[test]
    fn the_product_matches_a_central_difference_in_mach() {
        let step = 1.0e-6;
        for symmetry in symmetries() {
            let mach = 0.45;
            let solver = solver(mach, symmetry);
            let plus = self::solver(mach + step, symmetry);
            let minus = self::solver(mach - step, symmetry);
            let gamma = circulation(solver.mesh().vortex_count());
            for at_bound in [false, true] {
                let matrix = |s: &VortexLatticeSolver| {
                    if at_bound {
                        s.velocity_matrix().clone()
                    } else {
                        s.control_velocity_matrix().clone()
                    }
                };
                let finite = (matrix(&plus) - matrix(&minus)) * &gamma / (2.0 * step);
                let analytic = velocity_influence_mach_derivative_times(&solver, at_bound, &gamma);
                let scale = finite.amax();
                assert!(scale > 0.0);
                assert!(
                    (analytic - &finite).amax() <= 1.0e-6 * scale,
                    "symmetry {symmetry:?}, at_bound {at_bound}"
                );
            }
        }
    }

    #[test]
    fn incompressible_mach_skips_the_lattice() {
        for mach in [0.0, 1.0e-14] {
            let solver = solver(mach, SymmetryConfig::default());
            let n = solver.mesh().vortex_count();
            let gamma = circulation(n);
            for at_bound in [false, true] {
                let product = velocity_influence_mach_derivative_times(&solver, at_bound, &gamma);
                assert_eq!(product, VecX::zeros(3 * n));
                assert_eq!(
                    velocity_influence_mach_derivative(&solver, at_bound),
                    MatX::zeros(3 * n, n)
                );
            }
        }
    }
}
