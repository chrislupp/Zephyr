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

//! Rejection paths of the public API: each call below must refuse its input
//! with an error rather than return a partial or silently wrong result.

mod common;

use common::{controlled_wing, rectangular_mesh, rectangular_wing};
use zephyr_vlm::aero::AerodynamicResults;
use zephyr_vlm::analysis::{Analysis, AnalysisRequest, GradientRequest, Input, Output};
use zephyr_vlm::core::{SymmetryType, Vec3};
use zephyr_vlm::error::ZephyrError;
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{
    GeometryVariable, GradientCalculator, StripLoadAdjoint, StripLoadCotangent,
    trefftz_geometry_jvp,
};
use zephyr_vlm::io::{AvlReader, MassReader};
use zephyr_vlm::solver::{FlightCondition, SymmetryConfig, VortexLatticeSolver};
use zephyr_vlm::visualization::{VtkExportOptions, VtkExporter};

fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 0.05,
        ..FlightCondition::default()
    }
}

fn deflected() -> FlightCondition {
    FlightCondition {
        control_deflections: vec![("elevator".into(), 0.02)],
        ..condition()
    }
}

fn gradient(output: Output, input: Input) -> AnalysisRequest {
    AnalysisRequest {
        outputs: vec![output],
        gradients: vec![GradientRequest { output, input }],
        ..AnalysisRequest::default()
    }
}

fn solved(aircraft: &Aircraft, c: &FlightCondition) -> VortexLatticeSolver {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    solver.solve(c).unwrap();
    solver
}

fn polar_wing() -> Aircraft {
    let mut wing = Surface::new("Polar wing");
    wing.set_chordwise_panel_count(2);
    wing.set_spanwise_panel_count(4);
    for y in [-2.0, 2.0] {
        let mut section = Section::new(Vec3::new(0.0, y, 0.0), 1.0, 0.0);
        section.set_viscous_polar([0.0, 0.01, 0.5, 0.008, 1.0, 0.012]);
        wing.add_section(section);
    }
    with_surface(wing)
}

fn with_surface(surface: Surface) -> Aircraft {
    let mut aircraft = Aircraft::new("Fixture");
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);
    aircraft.add_surface(surface);
    aircraft
}

#[test]
fn analysis_rejects_invalid_derivative_requests() {
    let mesh = rectangular_mesh(2, 4);
    let aircraft = rectangular_wing(2, 4);
    let c = condition();

    let mut on_mesh = Analysis::from_mesh(&mesh);
    for input in [
        Input::MomentReference(3),
        Input::MeshNodes(5),
        Input::Section {
            surface: 0,
            section: 0,
            variable: GeometryVariable::Chord,
        },
        Input::Control("missing".into()),
    ] {
        assert!(
            on_mesh
                .evaluate(&c, &gradient(Output::Lift, input))
                .is_err()
        );
    }

    let mut controlled = Analysis::from_aircraft(&controlled_wing());
    let section = Input::Section {
        surface: 0,
        section: 0,
        variable: GeometryVariable::Chord,
    };
    assert!(
        controlled
            .evaluate(&deflected(), &gradient(Output::Lift, section))
            .is_err()
    );
    let mut wing = Analysis::from_aircraft(&aircraft);
    assert!(
        wing.evaluate(&c, &gradient(Output::Lift, Input::MeshNodes(0)))
            .is_err()
    );
}

#[test]
fn analysis_rejects_unknown_control_names_like_the_solver() {
    let mesh = rectangular_mesh(2, 4);
    let mut misspelled = condition();
    misspelled.control_deflections.push(("typo".into(), 0.3));
    let request = AnalysisRequest {
        outputs: vec![Output::Lift],
        ..AnalysisRequest::default()
    };

    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&mesh).unwrap();
    assert!(matches!(
        solver.solve(&misspelled),
        Err(ZephyrError::InvalidFlightCondition(_))
    ));
    let mut analysis = Analysis::from_mesh(&mesh);
    assert!(matches!(
        analysis.evaluate(&misspelled, &request),
        Err(ZephyrError::InvalidFlightCondition(_))
    ));
    // The rejection must not poison the context for a valid condition.
    assert!(analysis.evaluate(&condition(), &request).is_ok());
}

#[test]
fn analysis_rejects_degenerate_references_and_empty_lattices() {
    // Aircraft sources skip mesh-input validation, so the analysis checks its
    // own reference quantities and lattice.
    let request = AnalysisRequest::all_outputs();
    let mut zero_area = rectangular_wing(2, 4);
    zero_area.set_reference_area(0.0);
    assert!(
        Analysis::from_aircraft(&zero_area)
            .evaluate(&condition(), &request)
            .is_err()
    );
    let mut no_surfaces = Aircraft::new("Empty");
    no_surfaces.set_reference_area(1.0);
    no_surfaces.set_reference_chord(1.0);
    no_surfaces.set_reference_span(1.0);
    assert!(
        Analysis::from_aircraft(&no_surfaces)
            .evaluate(&condition(), &request)
            .is_err()
    );
}

#[test]
fn coincident_surfaces_are_reported_as_singular() {
    // Two copies of one surface give identical influence rows.
    let mut twice = rectangular_wing(2, 4);
    let copy = twice.surfaces()[0].clone();
    twice.add_surface(copy);
    let c = condition();

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&twice).unwrap();
    assert!(matches!(solver.solve(&c), Err(ZephyrError::SingularSystem)));
    assert!(matches!(
        Analysis::from_aircraft(&twice).evaluate(&c, &AnalysisRequest::all_outputs()),
        Err(ZephyrError::SingularSystem)
    ));
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&twice);
    assert!(matches!(
        calculator.compute_alpha_gradients(&c),
        Err(ZephyrError::SingularSystem)
    ));
}

#[test]
fn gradient_calculator_requires_a_matching_geometry_source() {
    let c = condition();
    let unset = GradientCalculator::new();
    assert!(unset.compute_alpha_gradients(&c).is_err());

    let mut on_mesh = GradientCalculator::new();
    on_mesh.set_mesh(&rectangular_mesh(2, 4)).unwrap();
    assert!(on_mesh.compute_geometry_gradients(&c).is_err());
    assert!(
        on_mesh
            .compute_section_gradient(&c, 0, 0, GeometryVariable::Chord)
            .is_err()
    );

    let mut on_aircraft = GradientCalculator::new();
    on_aircraft.set_aircraft(&rectangular_wing(2, 4));
    assert!(on_aircraft.compute_mesh_node_gradients(&c).is_err());
    assert!(
        on_aircraft
            .compute_section_gradient(&c, 0, 7, GeometryVariable::Chord)
            .is_err()
    );
}

#[test]
fn geometry_adjoints_reject_unsupported_solutions() {
    let mut on_mesh = GradientCalculator::new();
    let mesh = rectangular_mesh(2, 4);
    on_mesh.set_mesh(&mesh).unwrap();
    let mut unsolved = VortexLatticeSolver::new();
    unsolved.set_mesh(&mesh).unwrap();
    assert!(
        on_mesh
            .mesh_node_gradients_from(&unsolved, &condition())
            .is_err()
    );

    // Deflected control normals were once rejected here. The geometry adjoint
    // now differentiates through them (checked against re-solves in
    // tests/mesh_controls.rs), so the deflected solve must succeed.
    let mut controlled = GradientCalculator::new();
    controlled.set_aircraft(&controlled_wing());
    assert!(controlled.compute_geometry_gradients(&deflected()).is_ok());

    let mut symmetric = VortexLatticeSolver::new();
    symmetric.set_mesh(&mesh).unwrap();
    symmetric.set_symmetry(SymmetryConfig {
        z_symmetry: SymmetryType::Symmetric,
        z_symmetry_plane: -5.0,
        ..SymmetryConfig::default()
    });
    symmetric.solve(&condition()).unwrap();
    assert!(
        on_mesh
            .mesh_node_gradients_from(&symmetric, &condition())
            .is_err()
    );
}

#[test]
fn strip_load_adjoint_rejects_unsupported_solutions() {
    let points = [Vec3::zeros()];
    let mut unsolved = VortexLatticeSolver::new();
    unsolved.set_aircraft(&rectangular_wing(2, 4)).unwrap();
    assert!(StripLoadAdjoint::new(&unsolved, &condition(), &points).is_err());

    let deflected = deflected();
    let controlled = solved(&controlled_wing(), &deflected);
    assert!(StripLoadAdjoint::new(&controlled, &deflected, &points).is_err());

    let c = condition();
    let polar = solved(&polar_wing(), &c);
    assert!(StripLoadAdjoint::new(&polar, &c, &points).is_err());

    let mesh = rectangular_mesh(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&mesh).unwrap();
    solver.solve(&c).unwrap();
    let strips = solver.mesh().strip_count();
    let adjoint = StripLoadAdjoint::new(&solver, &c, &vec![Vec3::zeros(); strips]).unwrap();
    let short = StripLoadCotangent::zeros(strips - 1);
    assert!(adjoint.mesh_node_cotangents(&mesh, &[short]).is_err());
    let other = rectangular_mesh(3, 4);
    let column = StripLoadCotangent::zeros(strips);
    assert!(adjoint.mesh_node_cotangents(&other, &[column]).is_err());
}

#[test]
fn mesh_adjoints_reject_a_mirror_the_lattice_does_not_have() {
    // The solved mesh is not duplicated about y = 0; the mesh handed to the
    // adjoint claims it is, so the mirrored lattice it expects is absent.
    let c = condition();
    let mut solved_mesh = rectangular_mesh(2, 4);
    solved_mesh.surfaces[0].y_duplicated = false;
    let mut claimed = solved_mesh.clone();
    claimed.surfaces[0].y_duplicated = true;

    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&solved_mesh).unwrap();
    solver.solve(&c).unwrap();
    let strips = solver.mesh().strip_count();
    let adjoint = StripLoadAdjoint::new(&solver, &c, &vec![Vec3::zeros(); strips]).unwrap();
    let column = StripLoadCotangent::zeros(strips);
    assert!(adjoint.mesh_node_cotangents(&claimed, &[column]).is_err());

    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&claimed).unwrap();
    assert!(calculator.mesh_node_gradients_from(&solver, &c).is_err());
}

#[test]
fn trefftz_jvp_rejects_mismatched_directions() {
    let c = condition();
    let solver = solved(&rectangular_wing(2, 4), &c);
    let mesh = solver.mesh();
    let n = mesh.vortex_count();
    let result = trefftz_geometry_jvp(
        mesh,
        solver.solution(),
        4.0,
        0.0,
        0.0,
        &vec![Vec3::zeros(); n],
        &vec![Vec3::zeros(); n],
        &[],
    );
    assert!(result.is_err());
}

#[test]
fn stability_and_control_derivatives_require_a_valid_solution() {
    let aircraft = controlled_wing();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    assert!(
        solver
            .control_circulation_derivative(&condition(), "elevator")
            .is_err()
    );
}

#[test]
fn lattice_generation_rejects_degenerate_surfaces() {
    let mut one_section = Surface::new("Stub");
    one_section.set_chordwise_panel_count(2);
    one_section.set_spanwise_panel_count(2);
    one_section.add_section(Section::new(Vec3::zeros(), 1.0, 0.0));

    let mut no_chordwise = Surface::new("Flat");
    no_chordwise.set_chordwise_panel_count(0);
    no_chordwise.set_spanwise_panel_count(2);
    no_chordwise.add_section(Section::new(Vec3::new(0.0, -1.0, 0.0), 1.0, 0.0));
    no_chordwise.add_section(Section::new(Vec3::new(0.0, 1.0, 0.0), 1.0, 0.0));

    // Separated, but by less than the span tolerance.
    let mut negligible_span = Surface::new("Sliver");
    negligible_span.set_chordwise_panel_count(2);
    negligible_span.set_spanwise_panel_count(2);
    negligible_span.add_section(Section::new(Vec3::zeros(), 1.0, 0.0));
    negligible_span.add_section(Section::new(Vec3::new(0.0, 1e-14, 0.0), 1.0, 0.0));

    let mut zero_span = Surface::new("Point");
    zero_span.set_chordwise_panel_count(2);
    zero_span.set_spanwise_panel_count(2);
    zero_span.add_section(Section::new(Vec3::zeros(), 1.0, 0.0));
    zero_span.add_section(Section::new(Vec3::new(0.5, 0.0, 0.0), 1.0, 0.0));

    for surface in [one_section, no_chordwise, negligible_span, zero_span] {
        assert!(
            VortexLatticeSolver::new()
                .set_aircraft(&with_surface(surface))
                .is_err()
        );
    }
}

#[test]
fn readers_report_missing_files_and_truncated_input() {
    assert!(MassReader::read_path("does/not/exist.mass").is_err());
    let truncated = "Truncated\n0.0\n0 0 0\n1.0 1.0 1.0\n0.25 0.0 0.0\nSURFACE\n";
    assert!(AvlReader::default().read_str(truncated).is_err());
}

#[test]
fn vtk_exports_reject_unsolved_inputs_and_bad_paths() {
    let options = VtkExportOptions::default();
    let c = condition();
    let mut unsolved = VortexLatticeSolver::new();
    unsolved.set_aircraft(&rectangular_wing(2, 4)).unwrap();
    let dir = tempfile::tempdir().unwrap();

    assert!(
        VtkExporter::export_spanwise_loads(
            dir.path().join("loads.vtp"),
            unsolved.mesh(),
            &AerodynamicResults::default(),
            options,
        )
        .is_err()
    );
    assert!(
        VtkExporter::export_wake(
            dir.path().join("wake.vtp"),
            unsolved.mesh(),
            unsolved.solution(),
            &c,
            1.0,
            options,
        )
        .is_err()
    );
    assert!(
        VtkExporter::export_velocity_field(
            dir.path().join("field.vtp"),
            unsolved.mesh(),
            unsolved.solution(),
            &c,
            options,
        )
        .is_err()
    );

    let solver = solved(&rectangular_wing(2, 4), &c);
    for path in [dir.path().join(".."), dir.path().join("no_extension")] {
        assert!(
            VtkExporter::export_surface_mesh(&path, solver.mesh(), solver.solution(), options)
                .is_err()
        );
    }
}
