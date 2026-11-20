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

//! Ported from the frozen v0.2.0 `vtk_export_test.cpp`.
//!
//! Where the C++ returned `bool` and skipped when VTK was unavailable, the Rust
//! returns `Result` and gates on the `vtk` feature, so success is asserted
//! through `is_ok` and the feature-off behaviour lives in
//! `tests/visualization_disabled.rs`. Exports go to a `tempfile` directory
//! rather than the working directory, so the cases need no teardown and cannot
//! collide when run in parallel.

#![cfg(feature = "vtk")]

use zephyr_vlm::aero::{AerodynamicResults, ForceCalculator};
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, SolutionState, VortexLatticeSolver};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter, is_available};

/// The frozen fixture, solved at five degrees.
fn fixture() -> (
    VortexLatticeSolver,
    SolutionState,
    AerodynamicResults,
    FlightCondition,
    f64,
) {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(6);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Test Aircraft");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(2.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    let results = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    let chord = aircraft.reference_chord();
    (solver, solution, results, condition, chord)
}

/// Ports `IsAvailableReturnsConsistentValue`.
///
/// With the `vtk` feature on, availability is `true`; the feature-off case is
/// asserted in `tests/visualization_disabled.rs`.
#[test]
fn availability_matches_the_enabled_feature() {
    assert!(is_available());
}

/// Ports `ExportSurfaceMeshCreatesFile`.
#[test]
fn exporting_the_surface_mesh_writes_a_non_empty_file() {
    let (solver, solution, _, _, _) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("surface.vtu");
    VtkExporter::export_surface_mesh(&path, solver.mesh(), &solution, VtkExportOptions::default())
        .unwrap();
    assert!(path.exists());
    assert!(std::fs::metadata(&path).unwrap().len() > 0);
}

/// Ports `ExportSpanwiseDistributionsCreatesFile`.
#[test]
fn exporting_spanwise_loads_writes_a_file() {
    let (solver, _, results, _, _) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("loads.vtp");
    VtkExporter::export_spanwise_loads(&path, solver.mesh(), &results, VtkExportOptions::default())
        .unwrap();
    assert!(path.exists());
}

/// Ports `ExportWakeGeometryCreatesFile`.
#[test]
fn exporting_the_wake_writes_a_file() {
    let (solver, solution, _, condition, chord) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("wake.vtp");
    VtkExporter::export_wake(
        &path,
        solver.mesh(),
        &solution,
        &condition,
        chord,
        VtkExportOptions {
            wake_length_chords: 5.0,
            ..VtkExportOptions::default()
        },
    )
    .unwrap();
    assert!(path.exists());
}

/// Ports `ExportVelocityFieldCreatesFile`.
#[test]
fn exporting_the_velocity_field_writes_a_file() {
    let (solver, solution, _, condition, _) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("velocity.vtp");
    VtkExporter::export_velocity_field(
        &path,
        solver.mesh(),
        &solution,
        &condition,
        VtkExportOptions::default(),
    )
    .unwrap();
    assert!(path.exists());
}

/// Ports `ExportAllCreatesAllFiles`.
///
/// The four product filenames are part of the contract: a downstream ParaView
/// state file refers to them by name.
#[test]
fn exporting_everything_writes_all_four_products() {
    let (solver, solution, results, condition, chord) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("vtk_output");
    VtkExporter::export_all(
        &output,
        solver.mesh(),
        &solution,
        &results,
        &condition,
        chord,
        VtkExportOptions::default(),
    )
    .unwrap();
    for name in [
        "surface_mesh.vtu",
        "spanwise_loads.vtp",
        "wake.vtp",
        "velocity_field.vtp",
    ] {
        assert!(output.join(name).exists(), "{name} was not written");
    }
}

/// Ports `BinaryVsAsciiFormat`.
#[test]
fn binary_output_is_smaller_than_ascii() {
    let (solver, solution, _, _, _) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let mut sizes = Vec::new();
    for binary in [true, false] {
        let path = directory
            .path()
            .join(if binary { "b.vtu" } else { "a.vtu" });
        VtkExporter::export_surface_mesh(
            &path,
            solver.mesh(),
            &solution,
            VtkExportOptions {
                binary_format: binary,
                ..VtkExportOptions::default()
            },
        )
        .unwrap();
        sizes.push(std::fs::metadata(&path).unwrap().len());
    }
    assert!(
        sizes[0] < sizes[1],
        "binary {} >= ascii {}",
        sizes[0],
        sizes[1]
    );
}

/// Ports `InvalidSolutionReturnsFalse`.
///
/// An invalid solution must be refused *and* leave no partial file behind — a
/// truncated `.vtu` on disk is worse than none.
#[test]
fn an_invalid_solution_is_refused_and_writes_nothing() {
    let (solver, _, _, _, _) = fixture();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("invalid.vtu");
    let result = VtkExporter::export_surface_mesh(
        &path,
        solver.mesh(),
        &SolutionState::default(),
        VtkExportOptions::default(),
    );
    assert!(result.is_err());
    assert!(!path.exists());
}

/// Ports `ExportOptionsDefaultsAreReasonable`.
///
/// Pressure, circulation, wake, induced velocity and spanwise loads are on by
/// default; normals and indices are diagnostic and stay off. Binary is the
/// default because these files get large.
#[test]
fn export_option_defaults() {
    let options = VtkExportOptions::default();
    assert!(options.export_cp);
    assert!(options.export_circulation);
    assert!(!options.export_normals);
    assert!(!options.export_indices);
    assert!(options.export_wake);
    assert!(options.wake_length_chords > 0.0);
    assert!(options.export_induced_velocity);
    assert!(options.export_spanwise_distributions);
    assert!(options.binary_format);
}
