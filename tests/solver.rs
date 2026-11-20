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

mod common;

use approx::assert_abs_diff_eq;
use zephyr_vlm::ZephyrError;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, SpacingType, SymmetryType, Vec3, VecX};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{
    FlightCondition, InfluenceMatrixBuilder, SymmetryConfig, VortexLattice, VortexLatticeSolver,
    horseshoe_vortex_velocity,
};

use common::{rectangular_mesh, rectangular_wing};

#[test]
fn aircraft_lattice_preserves_expected_panel_order_and_counts() {
    let aircraft = rectangular_wing(3, 5);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert_eq!(lattice.vortex_count(), 15);
    assert_eq!(lattice.strip_count(), 5);
    assert_eq!(lattice.surface_count(), 1);
    for (strip_index, strip) in lattice.strips().iter().enumerate() {
        assert_eq!(strip.first_vortex_index, strip_index * 3);
        assert_eq!(strip.num_vortices, 3);
        for offset in 0..3 {
            assert_eq!(
                lattice.vortices()[strip.first_vortex_index + offset].strip_index,
                strip_index
            );
        }
    }
}

#[test]
fn flat_rectangular_lattice_has_quarter_chord_vortices_and_three_quarter_controls() {
    let aircraft = rectangular_wing(1, 2);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    for vortex in lattice.vortices() {
        assert_abs_diff_eq!(vortex.bound_vortex_center.x, 0.25, epsilon = 1.0e-12);
        assert_abs_diff_eq!(vortex.control_point.x, 0.75, epsilon = 1.0e-12);
        assert_abs_diff_eq!(vortex.normal_at_control_point, Vec3::z(), epsilon = 1.0e-12);
    }
}

#[test]
fn twist_rotates_panel_normals() {
    let mut wing = Surface::new("Twisted");
    wing.set_chordwise_panel_count(1);
    wing.set_spanwise_panel_count(2);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 10.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 10.0));
    let mut aircraft = Aircraft::new("Twisted");
    aircraft.add_surface(wing);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert!(lattice.vortices()[0].normal_at_control_point.x > 0.1);
    assert_abs_diff_eq!(
        lattice.vortices()[0].normal_at_control_point.norm(),
        1.0,
        epsilon = 1.0e-12
    );
}

#[test]
fn y_duplicate_mirrors_geometry_and_preserves_orientation() {
    let mut aircraft = rectangular_wing(2, 3);
    let wing = &mut aircraft.surfaces_mut()[0];
    wing.set_y_duplicate(0.0);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    assert_eq!(lattice.surface_count(), 2);
    assert_eq!(lattice.vortex_count(), 12);
    let original = &lattice.surface_info()[0];
    let mirror = &lattice.surface_info()[1];
    assert!(mirror.is_y_duplicated);
    for offset in 0..original.num_vortices {
        let left = &lattice.vortices()[original.first_vortex_index + offset];
        let right = &lattice.vortices()[mirror.first_vortex_index + offset];
        assert_abs_diff_eq!(left.control_point.x, right.control_point.x);
        assert_abs_diff_eq!(left.control_point.y, -right.control_point.y);
        assert!(right.bound_vortex_right.y > right.bound_vortex_left.y);
    }
}

#[test]
fn direct_mesh_and_aircraft_paths_match_for_rectangular_geometry() {
    let mut aircraft = rectangular_wing(2, 4);
    aircraft.surfaces_mut()[0].set_chordwise_spacing(SpacingType::Uniform);
    aircraft.surfaces_mut()[0].set_spanwise_spacing(SpacingType::Uniform);
    let mesh = rectangular_mesh(2, 4);
    let mut aircraft_solver = VortexLatticeSolver::new();
    let mut mesh_solver = VortexLatticeSolver::new();
    aircraft_solver.set_aircraft(&aircraft).unwrap();
    mesh_solver.set_mesh(&mesh).unwrap();
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let aircraft_solution = aircraft_solver.solve(&condition).unwrap();
    let mesh_solution = mesh_solver.solve(&condition).unwrap();
    assert_abs_diff_eq!(
        aircraft_solution.circulation,
        mesh_solution.circulation,
        epsilon = 1.0e-11
    );
}

#[test]
fn direct_mesh_dimensions_and_normal_recalculation_are_valid() {
    let mesh = rectangular_mesh(2, 3);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_mesh(&mesh).unwrap();
    assert_eq!(lattice.vortex_count(), 6);
    assert_eq!(lattice.strip_count(), 3);
    assert!(lattice.mesh_input().is_some());
    lattice.recalculate_normals();
    assert!(
        lattice
            .vortices()
            .iter()
            .all(|vortex| (vortex.normal_at_control_point.norm() - 1.0).abs() < 1.0e-12)
    );
}

#[test]
fn regularized_horseshoe_kernel_is_finite_near_the_vortex() {
    let left = Vec3::new(0.0, -0.5, 0.0);
    let right = Vec3::new(0.0, 0.5, 0.0);
    for point in [
        Vec3::new(0.0, 0.0, 1.0e-10),
        Vec3::new(1.0, 0.0, 0.2),
        Vec3::new(-1.0, 0.2, 0.1),
    ] {
        let velocity = horseshoe_vortex_velocity(&point, &left, &right, 0.9, 0.01, true);
        assert!(velocity.iter().all(|value| value.is_finite()));
    }
}

#[test]
fn excluding_bound_segment_changes_only_the_finite_segment_contribution() {
    let left = Vec3::new(0.0, -0.5, 0.0);
    let right = Vec3::new(0.0, 0.5, 0.0);
    let point = Vec3::new(0.5, 0.0, 0.2);
    let with_bound = horseshoe_vortex_velocity(&point, &left, &right, 1.0, 0.001, true);
    let without_bound = horseshoe_vortex_velocity(&point, &left, &right, 1.0, 0.001, false);
    assert!((with_bound - without_bound).norm() > 1.0e-4);
}

#[test]
fn influence_matrix_dimensions_and_normal_projection_are_consistent() {
    let aircraft = rectangular_wing(2, 4);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();
    let builder = InfluenceMatrixBuilder::new(&lattice);
    let (aic, velocity) = builder.normalwash_and_velocity();
    assert_eq!(aic.shape(), (8, 8));
    assert_eq!(velocity.shape(), (24, 8));
    for row in 0..8 {
        for column in 0..8 {
            let normal = lattice.vortices()[row].normal_at_control_point;
            let projected = normal.dot(&builder.single_vortex_velocity(
                &lattice.vortices()[row].control_point,
                column,
                true,
            ));
            assert_abs_diff_eq!(aic[(row, column)], projected, epsilon = 1.0e-14);
        }
    }
}

#[test]
fn aic_is_nonsingular_and_unit_response_residual_is_small() {
    let aircraft = rectangular_wing(3, 6);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve_unit_responses().unwrap();
    assert!(solver.is_aic_valid());
    assert!(solver.are_unit_responses_valid());
    let solution = solver.solve(&FlightCondition::default()).unwrap();
    let rhs = VecX::from_iterator(
        solver.mesh().vortex_count(),
        solver
            .mesh()
            .vortices()
            .iter()
            .map(|vortex| -vortex.normal_at_control_point.z),
    );
    let residual = solver.aic() * &solution.circulation_per_w - rhs;
    assert!(residual.norm() < 1.0e-10);
}

#[test]
fn zero_alpha_flat_wing_has_zero_circulation() {
    let aircraft = rectangular_wing(3, 6);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&FlightCondition::default()).unwrap();
    assert!(solution.circulation.norm() < 1.0e-13);
    assert!(solution.induced_velocity.norm() < 1.0e-13);
}

#[test]
fn positive_and_negative_alpha_reverse_circulation() {
    let aircraft = rectangular_wing(3, 6);
    let solve = |alpha| {
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft).unwrap();
        solver
            .solve(&FlightCondition {
                alpha,
                ..FlightCondition::default()
            })
            .unwrap()
            .circulation
    };
    let positive = solve(3.0 * DEGREES_TO_RADIANS);
    let negative = solve(-3.0 * DEGREES_TO_RADIANS);
    assert_abs_diff_eq!(positive, -negative, epsilon = 1.0e-12);
}

#[test]
fn freestream_direction_handles_alpha_and_beta() {
    let direction = VortexLatticeSolver::freestream_direction(&FlightCondition {
        alpha: 10.0 * DEGREES_TO_RADIANS,
        beta: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    });
    assert_abs_diff_eq!(direction.norm(), 1.0, epsilon = 1.0e-14);
    assert!(direction.x > 0.9);
    assert!(direction.y > 0.0);
    assert!(direction.z > 0.0);
}

#[test]
fn mach_change_rebuilds_the_cached_aic() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    let incompressible = solver.aic().clone();
    solver
        .solve(&FlightCondition {
            mach: 0.5,
            ..FlightCondition::default()
        })
        .unwrap();
    assert_eq!(solver.current_mach(), 0.5);
    assert!((&incompressible - solver.aic()).norm() > 1.0e-5);
}

/// Transonic and supersonic conditions are rejected rather than answered.
///
/// The builder clamps the Prandtl-Glauert factor to 1 at and above Mach 1, so
/// `solve` used to return the *incompressible* solution for `mach = 1.5` — a
/// finite, plausible-looking answer to a question the method cannot address,
/// with `mach_number()` still reporting 1.5.
#[test]
fn supersonic_mach_is_rejected() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    for mach in [1.0, 1.5, 3.0] {
        let error = solver
            .solve(&FlightCondition {
                mach,
                ..FlightCondition::default()
            })
            .unwrap_err();
        assert!(
            matches!(error, ZephyrError::InvalidFlightCondition(_)),
            "expected Mach {mach} to be rejected, got {error:?}"
        );
    }
    // Just below Mach 1 is still a valid compressible solve.
    let solution = solver
        .solve(&FlightCondition {
            mach: 0.95,
            ..FlightCondition::default()
        })
        .unwrap();
    assert!(solution.is_valid);
}

/// Non-finite body rates must be rejected, not solved.
///
/// Only `alpha`, `beta`, `mach` and the control deflections were checked, so a
/// NaN roll rate produced an all-NaN circulation reported as valid. NaN also
/// defeats `dynamic_pressure < 0.0`, since every ordered comparison against NaN
/// is false.
#[test]
fn non_finite_rates_and_dynamic_pressure_are_rejected() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();

    let bad = [
        FlightCondition {
            roll_rate: f64::NAN,
            ..FlightCondition::default()
        },
        FlightCondition {
            pitch_rate: f64::INFINITY,
            ..FlightCondition::default()
        },
        FlightCondition {
            yaw_rate: f64::NAN,
            ..FlightCondition::default()
        },
        FlightCondition {
            dynamic_pressure: f64::NAN,
            ..FlightCondition::default()
        },
    ];
    for condition in bad {
        let error = solver.solve(&condition).unwrap_err();
        assert!(
            matches!(error, ZephyrError::InvalidFlightCondition(_)),
            "expected rejection, got {error:?}"
        );
    }
}

#[test]
fn symmetry_images_change_the_influence_matrix() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.rebuild_aic().unwrap();
    let base = solver.aic().clone();
    solver.set_symmetry(SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        y_symmetry_plane: 0.0,
        ..SymmetryConfig::default()
    });
    solver.rebuild_aic().unwrap();
    assert!((&base - solver.aic()).norm() > 1.0e-6);
}

#[test]
fn mesh_mutation_invalidates_dependent_caches() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    assert!(solver.is_aic_valid());
    solver.mutate_mesh(|mesh| {
        // Any mesh mutation must drop the cached factorization. Regenerating
        // is enough and, unlike the core-radius knobs this test used to poke,
        // genuinely rebuilds the system.
        mesh.generate_from_aircraft(&rectangular_wing(2, 4))
            .unwrap();
    });
    assert!(!solver.is_aic_valid());
    assert!(!solver.are_unit_responses_valid());
    assert!(!solver.solution().is_valid);
}

/// The core-radius fractions that actually reach the kernel are the builder's.
///
/// `VortexLattice` used to carry a second, independent pair of fractions with
/// different defaults (0.75/0.25 against the builder's 0.0001/0.0001) that were
/// never passed to assembly, so setting them left the AIC bit-identical while
/// appearing to tune the regularization. They have been removed; this pins the
/// remaining knobs to being live.
///
/// The fractions only apply between panels on *different* components — a
/// same-component pair uses the small fixed `1e-4 * panel_width` core — so the
/// fixture needs two surfaces with distinct component indices. A single-surface
/// aircraft would leave the branch unreachable and the test vacuous.
#[test]
fn builder_core_radius_fractions_change_the_assembled_matrix() {
    let mut wing = Surface::new("wing");
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    wing.set_chordwise_panel_count(2);
    wing.set_spanwise_panel_count(3);
    wing.set_component_index(1);

    let mut tail = Surface::new("tail");
    tail.add_section(Section::new(Vec3::new(3.0, 0.0, 0.3), 0.6, 0.0));
    tail.add_section(Section::new(Vec3::new(3.0, 1.0, 0.3), 0.6, 0.0));
    tail.set_chordwise_panel_count(2);
    tail.set_spanwise_panel_count(2);
    tail.set_component_index(2);

    let mut aircraft = Aircraft::new("two components");
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();

    let assemble = |chord_fraction: f64, width_fraction: f64| -> Vec<u64> {
        let mut builder = InfluenceMatrixBuilder::new(solver.mesh());
        builder.set_core_radius_chord_fraction(chord_fraction);
        builder.set_core_radius_width_fraction(width_fraction);
        builder
            .normalwash_matrix()
            .iter()
            .map(|value| value.to_bits())
            .collect()
    };

    let baseline = assemble(0.0001, 0.0001);
    assert_ne!(baseline, assemble(0.3, 0.0001));
    assert_ne!(baseline, assemble(0.0001, 0.3));
    assert_eq!(baseline, assemble(0.0001, 0.0001));
}

#[test]
fn invalid_flight_conditions_return_errors() {
    let aircraft = rectangular_wing(1, 2);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    for condition in [
        FlightCondition {
            alpha: f64::NAN,
            ..FlightCondition::default()
        },
        FlightCondition {
            mach: -0.1,
            ..FlightCondition::default()
        },
        FlightCondition {
            dynamic_pressure: -1.0,
            ..FlightCondition::default()
        },
        FlightCondition {
            control_deflections: vec![("bad".into(), f64::INFINITY)],
            ..FlightCondition::default()
        },
    ] {
        assert!(solver.solve(&condition).is_err());
    }
}

#[test]
fn empty_solver_rejects_solve_and_aic_build() {
    let mut solver = VortexLatticeSolver::new();
    assert!(solver.rebuild_aic().is_err());
    assert!(solver.solve(&FlightCondition::default()).is_err());
}

/// The influence kernel must not depend on the unit the model is expressed in.
///
/// The bound-leg and trailing-leg guards used to compare a length⁴ and a
/// length quantity against an absolute `EPSILON`, so below roughly millimetre
/// scale the bound contribution was silently dropped from every panel. That
/// moved `‖AIC‖` by 3.3% at 1e-3 and 4.2% at 1e-4 with no error returned.
///
/// The AIC carries units of inverse length, so `scale · ‖AIC‖` is the
/// dimensionless invariant: it must be identical for a geometrically similar
/// wing at any scale.
#[test]
fn influence_assembly_is_invariant_under_uniform_scaling() {
    fn wing_at_scale(scale: f64) -> Aircraft {
        let mut surface = Surface::new("wing");
        surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), scale, 0.0));
        surface.add_section(Section::new(Vec3::new(0.0, 4.0 * scale, 0.0), scale, 0.0));
        surface.set_spanwise_panel_count(4);
        surface.set_chordwise_panel_count(2);
        let mut aircraft = Aircraft::new("scaled");
        aircraft.add_surface(surface);
        aircraft
    }

    fn invariant(scale: f64) -> f64 {
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&wing_at_scale(scale)).unwrap();
        solver.rebuild_aic().unwrap();
        scale * solver.aic().iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    let reference = invariant(1.0);
    assert!(reference.is_finite() && reference > 0.0);
    for scale in [1.0e3, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5] {
        assert_abs_diff_eq!(invariant(scale), reference, epsilon = 1.0e-9);
    }
}

/// `induced_velocity` must agree with the assembled influence matrix.
///
/// Assembly adds four images when both symmetry planes are active — y, z, and
/// the doubly-reflected yz — but the field-point helper omitted the last one.
/// With both planes on, the two disagreed by roughly 3% including a sign flip
/// in the y component, so a caller sampling the velocity field got different
/// physics from the one solving the system.
///
/// The core radii are pinned so the comparison is exact: assembly uses the
/// small same-component core, while `induced_velocity` always takes the
/// cross-component branch. Setting the chord fraction to zero and the width
/// fraction to the same 1e-4 makes both reduce to `1e-4 * panel_width`.
#[test]
fn induced_velocity_matches_the_influence_matrix_under_double_symmetry() {
    let aircraft = rectangular_wing(2, 3);
    let mut lattice = VortexLattice::default();
    lattice.generate_from_aircraft(&aircraft).unwrap();

    for symmetry in [
        SymmetryConfig {
            y_symmetry: SymmetryType::Symmetric,
            y_symmetry_plane: 0.0,
            z_symmetry: SymmetryType::Symmetric,
            z_symmetry_plane: 0.4,
        },
        SymmetryConfig {
            y_symmetry: SymmetryType::Symmetric,
            y_symmetry_plane: 0.0,
            z_symmetry: SymmetryType::Antisymmetric,
            z_symmetry_plane: 0.4,
        },
        SymmetryConfig {
            y_symmetry: SymmetryType::Antisymmetric,
            y_symmetry_plane: 0.0,
            z_symmetry: SymmetryType::Symmetric,
            z_symmetry_plane: 0.4,
        },
    ] {
        let mut builder = InfluenceMatrixBuilder::new(&lattice);
        builder.set_symmetry(symmetry);
        builder.set_core_radius_chord_fraction(0.0);
        builder.set_core_radius_width_fraction(0.0001);
        let velocity = builder.velocity_influence_matrix();

        let n = lattice.vortex_count();
        for source in 0..n {
            // Skip the diagonal: assembly drops the self-influence bound leg.
            let target = (source + 1) % n;
            let mut circulation = VecX::zeros(n);
            circulation[source] = 1.0;
            let sampled = builder.induced_velocity(
                &lattice.vortices()[target].bound_vortex_center,
                &circulation,
            );
            let assembled = Vec3::new(
                velocity[(3 * target, source)],
                velocity[(3 * target + 1, source)],
                velocity[(3 * target + 2, source)],
            );
            assert_abs_diff_eq!(sampled, assembled, epsilon = 1.0e-12);
        }
    }
}
