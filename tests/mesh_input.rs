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

//! Ported from the frozen v0.2.0 `mesh_input_test.cpp`.
//!
//! `MeshInput` is the direct-mesh entry point: node coordinates handed in from
//! outside rather than generated from sections. The point of these cases is
//! that it reaches the same answer as the `Aircraft` path for geometry the two
//! describe identically.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, SpacingType, Vec3};
use zephyr_vlm::geometry::{Aircraft, MeshInput, Section, Surface, SurfaceMesh};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn rectangular_mesh(nx: usize, ny: usize, chord: f64, span: f64) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", nx, ny);
    for i in 0..=nx {
        for j in 0..=ny {
            let x = chord * i as f64 / nx as f64;
            let y = span * j as f64 / ny as f64;
            surface.set_node(i, j, Vec3::new(x, y, 0.0)).unwrap();
        }
    }
    let mut input = MeshInput::default();
    input.surfaces.push(surface);
    input.reference.area = chord * span;
    input.reference.chord = chord;
    input.reference.span = span;
    input.reference.moment_reference_point = Vec3::new(0.25 * chord, 0.5 * span, 0.0);
    input
}

fn rectangular_aircraft(nx: usize, ny: usize, chord: f64, span: f64) -> Aircraft {
    let mut surface = Surface::new("Wing");
    surface.set_chordwise_panel_count(nx);
    surface.set_chordwise_spacing(SpacingType::Uniform);
    surface.set_spanwise_panel_count(ny);
    surface.set_spanwise_spacing(SpacingType::Uniform);
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), chord, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, span, 0.0), chord, 0.0));

    let mut aircraft = Aircraft::new("RectWing");
    aircraft.set_reference_area(chord * span);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft.set_moment_reference_point(Vec3::new(0.25 * chord, 0.5 * span, 0.0));
    aircraft.add_surface(surface);
    aircraft
}

/// Ports `BasicSolveProducesValidResults`.
///
/// The lift bounds are loose on purpose: they pin the result to the physically
/// plausible band for a AR-4 flat plate at 5 degrees rather than to a specific
/// number, which is what the frozen case asserted.
#[test]
fn direct_mesh_solve_produces_physically_plausible_lift() {
    let input = rectangular_mesh(4, 8, 1.0, 4.0);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    assert!(solution.is_valid);
    assert!(!solution.circulation.is_empty());

    let results = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        &condition,
    );
    assert!(results.is_valid);
    assert!(results.coefficients.lift > 0.2);
    assert!(results.coefficients.lift < 1.0);
}

/// Ports `ConsistencyWithAircraftPath`.
#[test]
fn direct_mesh_and_aircraft_paths_agree() {
    let (nx, ny, chord, span) = (6, 12, 1.0, 6.0);
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };

    let input = rectangular_mesh(nx, ny, chord, span);
    let mut mesh_solver = VortexLatticeSolver::new();
    mesh_solver.set_mesh(&input).unwrap();
    let mesh_solution = mesh_solver.solve(&condition).unwrap();
    let mesh_results = ForceCalculator::new().compute_with_references(
        &input.reference,
        mesh_solver.mesh(),
        &mesh_solution,
        &condition,
    );

    let aircraft = rectangular_aircraft(nx, ny, chord, span);
    let mut aircraft_solver = VortexLatticeSolver::new();
    aircraft_solver.set_aircraft(&aircraft).unwrap();
    let aircraft_solution = aircraft_solver.solve(&condition).unwrap();
    let aircraft_results = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        aircraft_solver.mesh(),
        &aircraft_solution,
        &condition,
    );

    assert!(mesh_results.is_valid);
    assert!(aircraft_results.is_valid);
    assert_abs_diff_eq!(
        mesh_results.coefficients.lift,
        aircraft_results.coefficients.lift,
        epsilon = 0.05
    );
    assert_abs_diff_eq!(
        mesh_results.coefficients.pitch,
        aircraft_results.coefficients.pitch,
        epsilon = 0.05
    );
}

/// Ports `MeshDimensions`.
#[test]
fn direct_mesh_lattice_dimensions_match_the_input() {
    let (nx, ny) = (3, 5);
    let input = rectangular_mesh(nx, ny, 1.0, 4.0);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&input).unwrap();
    let mesh = solver.mesh();
    assert_eq!(mesh.vortex_count(), nx * ny);
    assert_eq!(mesh.strip_count(), ny);
    assert_eq!(mesh.surface_info().len(), 1);
    assert_eq!(mesh.surface_info()[0].num_chordwise, nx);
    assert_eq!(mesh.surface_info()[0].num_spanwise, ny);
}
