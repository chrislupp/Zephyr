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

#![cfg(feature = "vtk")]

use approx::assert_abs_diff_eq;
use vtkio::model::{DataSet, IOBuffer, Piece};
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter};

fn solved_wing() -> (
    Aircraft,
    VortexLatticeSolver,
    zephyr_vlm::solver::SolutionState,
    zephyr_vlm::aero::AerodynamicResults,
    FlightCondition,
) {
    let mut surface = Surface::new("Wing");
    surface.set_chordwise_panel_count(2);
    surface.set_spanwise_panel_count(4);
    surface.add_section(Section::new(Vec3::new(0.0, -1.0, 0.0), 1.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("VTK fixture");
    aircraft.set_reference_area(2.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(2.0);
    aircraft.add_surface(surface);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: 0.08,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    let results = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    (aircraft, solver, solution, results, condition)
}

#[test]
fn all_four_default_products_round_trip_through_vtkio() {
    let (aircraft, solver, solution, results, condition) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    VtkExporter::export_all(
        directory.path(),
        solver.mesh(),
        &solution,
        &results,
        &condition,
        aircraft.reference_chord(),
        VtkExportOptions::default(),
    )
    .unwrap();
    for name in [
        "surface_mesh.vtu",
        "spanwise_loads.vtp",
        "wake.vtp",
        "velocity_field.vtp",
    ] {
        vtkio::Vtk::import(directory.path().join(name)).unwrap();
    }
}

#[test]
fn surface_arrays_and_topology_are_present() {
    let (_, solver, solution, _, _) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("surface.vtu");
    let options = VtkExportOptions {
        export_normals: true,
        export_indices: true,
        binary_format: false,
        ..VtkExportOptions::default()
    };
    VtkExporter::export_surface_mesh(&path, solver.mesh(), &solution, options).unwrap();
    let vtk = vtkio::Vtk::import(&path).unwrap();
    let DataSet::UnstructuredGrid { pieces, .. } = vtk.data else {
        panic!("expected an unstructured grid");
    };
    let Piece::Inline(piece) = &pieces[0] else {
        panic!("expected inline VTK data");
    };
    assert_eq!(piece.cells.num_cells(), solver.mesh().vortex_count());
    let names: Vec<_> = piece.data.cell.iter().map(|array| array.name()).collect();
    assert_eq!(names, ["Cp", "Circulation", "Surface", "Strip", "Normal"]);
}

#[test]
fn invalid_solution_creates_no_file() {
    let (_, solver, _, _, _) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("invalid.vtu");
    let result = VtkExporter::export_surface_mesh(
        &path,
        solver.mesh(),
        &zephyr_vlm::solver::SolutionState::default(),
        VtkExportOptions::default(),
    );
    assert!(result.is_err());
    assert!(!path.exists());
}

#[test]
fn binary_is_smaller_than_ascii_for_fixture() {
    let (_, solver, solution, _, _) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let binary = directory.path().join("binary.vtu");
    let ascii = directory.path().join("ascii.vtu");
    VtkExporter::export_surface_mesh(
        &binary,
        solver.mesh(),
        &solution,
        VtkExportOptions::default(),
    )
    .unwrap();
    VtkExporter::export_surface_mesh(
        &ascii,
        solver.mesh(),
        &solution,
        VtkExportOptions {
            binary_format: false,
            ..VtkExportOptions::default()
        },
    )
    .unwrap();
    assert!(std::fs::metadata(binary).unwrap().len() < std::fs::metadata(ascii).unwrap().len());
}

#[test]
fn vtk_availability_and_default_options_are_consistent() {
    assert!(zephyr_vlm::visualization::is_available());
    let options = VtkExportOptions::default();
    assert!(options.export_cp);
    assert!(options.export_circulation);
    assert!(options.export_wake);
    assert!(options.export_induced_velocity);
    assert!(options.export_spanwise_distributions);
    assert!(options.binary_format);
    assert!(options.wake_length_chords > 0.0);
}

#[test]
fn disabled_optional_products_are_omitted_by_export_all() {
    let (aircraft, solver, solution, results, condition) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    VtkExporter::export_all(
        directory.path(),
        solver.mesh(),
        &solution,
        &results,
        &condition,
        aircraft.reference_chord(),
        VtkExportOptions {
            export_wake: false,
            export_induced_velocity: false,
            export_spanwise_distributions: false,
            ..VtkExportOptions::default()
        },
    )
    .unwrap();
    assert!(directory.path().join("surface_mesh.vtu").exists());
    assert!(!directory.path().join("spanwise_loads.vtp").exists());
    assert!(!directory.path().join("wake.vtp").exists());
    assert!(!directory.path().join("velocity_field.vtp").exists());
}

#[test]
fn surface_output_omits_unselected_cell_arrays() {
    let (_, solver, solution, _, _) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("minimal.vtu");
    VtkExporter::export_surface_mesh(
        &path,
        solver.mesh(),
        &solution,
        VtkExportOptions {
            export_cp: false,
            export_circulation: false,
            export_normals: false,
            export_indices: false,
            binary_format: false,
            ..VtkExportOptions::default()
        },
    )
    .unwrap();
    let contents = std::fs::read_to_string(path).unwrap();
    for omitted in ["Cp", "Circulation", "Normal", "Surface", "Strip"] {
        assert!(!contents.contains(&format!("Name=\"{omitted}\"")));
    }
}

#[test]
fn wake_length_and_conditional_arrays_follow_options() {
    let (_, solver, solution, _, condition) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("wake.vtp");
    VtkExporter::export_wake(
        &path,
        solver.mesh(),
        &solution,
        &condition,
        2.0,
        VtkExportOptions {
            export_circulation: false,
            export_indices: false,
            wake_length_chords: 3.0,
            ..VtkExportOptions::default()
        },
    )
    .unwrap();
    let vtk = vtkio::Vtk::import(path).unwrap();
    let DataSet::PolyData { pieces, .. } = vtk.data else {
        panic!("expected polydata");
    };
    let Piece::Inline(piece) = &pieces[0] else {
        panic!("expected inline VTK data");
    };
    let IOBuffer::F64(points) = &piece.points else {
        panic!("expected f64 points");
    };
    let displacement = Vec3::new(
        points[3] - points[0],
        points[4] - points[1],
        points[5] - points[2],
    );
    assert_abs_diff_eq!(displacement.norm(), 6.0, epsilon = 1.0e-12);
    assert!(piece.data.cell.is_empty());
}

#[test]
fn spanwise_and_velocity_outputs_respect_array_flags() {
    let (_, solver, solution, results, condition) = solved_wing();
    let directory = tempfile::tempdir().unwrap();
    let span = directory.path().join("span.vtp");
    let velocity = directory.path().join("velocity.vtp");
    let options = VtkExportOptions {
        export_induced_velocity: false,
        export_indices: false,
        binary_format: false,
        ..VtkExportOptions::default()
    };
    VtkExporter::export_spanwise_loads(&span, solver.mesh(), &results, options).unwrap();
    VtkExporter::export_velocity_field(&velocity, solver.mesh(), &solution, &condition, options)
        .unwrap();
    let span_text = std::fs::read_to_string(span).unwrap();
    let velocity_text = std::fs::read_to_string(velocity).unwrap();
    assert!(!span_text.contains("Name=\"InducedVelocity\""));
    assert!(!span_text.contains("Name=\"Surface\""));
    assert!(!velocity_text.contains("Name=\"InducedVelocity\""));
    assert!(!velocity_text.contains("Name=\"Surface\""));
    assert!(velocity_text.contains("Name=\"TotalVelocity\""));
}
