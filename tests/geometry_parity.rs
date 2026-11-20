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

//! Ported from the frozen v0.2.0 `geometry_test.cpp`.
//!
//! The largest single file in the C++ suite: airfoil generation, section and
//! surface geometry, aircraft aggregation, and the AVL reader's surface and
//! section keywords. Where the C++ threw, the Rust returns a typed `Result`,
//! and those cases assert the error rather than the exception.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, SpacingType, SymmetryType, Vec3};
use zephyr_vlm::geometry::{
    Aircraft, Airfoil, Body, ControlSurfaceDefinition, DesignVariableDefinition, Section, Surface,
};
use zephyr_vlm::io::AvlReader;

const TOL: f64 = 1.0e-6;
/// The C++ built its NACA sections on a 61-point distribution.
const POINTS: usize = 61;

fn aileron() -> ControlSurfaceDefinition {
    ControlSurfaceDefinition {
        name: "aileron".into(),
        hinge_x_fraction: 0.75,
        hinge_vector: Vec3::new(0.0, 1.0, 0.0),
        gain: 1.0,
        sign: 1.0,
    }
}

// ============================================================================
// Airfoil
// ============================================================================

/// Ports `AirfoilTest.FlatPlateCreation`.
#[test]
fn flat_plate_creation() {
    let airfoil = Airfoil::flat_plate();
    assert!(airfoil.is_valid());
    assert_eq!(airfoil.name(), "Flat Plate");
    assert_abs_diff_eq!(airfoil.max_camber(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.0, epsilon = TOL);
    for x in [0.0, 0.5, 1.0] {
        assert_abs_diff_eq!(airfoil.camber(x), 0.0, epsilon = TOL);
    }
    assert_abs_diff_eq!(airfoil.camber_slope(0.5), 0.0, epsilon = TOL);
}

/// Ports `AirfoilTest.Naca0012Creation`.
///
/// The name is the code as written: `12`, not `0012`. A symmetric section has
/// zero camber everywhere, which is what separates it from 2412 below.
#[test]
fn naca_0012_creation() {
    let airfoil = Airfoil::naca4(12, POINTS).unwrap();
    assert!(airfoil.is_valid());
    assert_eq!(airfoil.name(), "NACA 12");
    assert_abs_diff_eq!(airfoil.max_camber(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.12, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.camber(0.3), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.camber(0.5), 0.0, epsilon = TOL);
}

/// Ports `AirfoilTest.Naca2412Creation`.
#[test]
fn naca_2412_creation() {
    let airfoil = Airfoil::naca4(2412, POINTS).unwrap();
    assert!(airfoil.is_valid());
    assert_eq!(airfoil.name(), "NACA 2412");
    assert_abs_diff_eq!(airfoil.max_camber(), 0.02, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_camber_position(), 0.4, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.12, epsilon = TOL);
    // The camber peak sits at 40% chord, which is the `4` in 2412.
    let peak = airfoil.camber(0.4);
    assert!(peak > 0.0);
    assert_abs_diff_eq!(peak, 0.02, epsilon = 0.001);
}

/// Ports `AirfoilTest.Naca4DigitFromString`.
#[test]
fn naca_4412_from_string() {
    let airfoil = Airfoil::naca4_str("4412", POINTS).unwrap();
    assert_eq!(airfoil.name(), "NACA 4412");
    assert_abs_diff_eq!(airfoil.max_camber(), 0.04, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_camber_position(), 0.4, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.12, epsilon = TOL);
}

/// Ports `AirfoilTest.CamberSlopeAtLeadingEdge`.
///
/// Camber rises from the leading edge to the 40% peak and falls after it, so
/// the slope changes sign across the chord. A flat or monotone camber line
/// would pass a positive-slope check alone.
#[test]
fn camber_slope_changes_sign_across_the_chord() {
    let airfoil = Airfoil::naca4(2412, POINTS).unwrap();
    assert!(airfoil.camber_slope(0.0) > 0.0);
    assert!(airfoil.camber_slope(1.0) < 0.0);
}

/// Ports `AirfoilExtraTest.InvalidNacaString`.
#[test]
fn invalid_naca_string_is_rejected() {
    assert!(Airfoil::naca4_str("xyz", POINTS).is_err());
}

/// Ports `AirfoilExtraTest.NacaCodeOutOfRange`.
///
/// The C++ checked for a negative code at runtime; `u16` rules that out at the
/// type level, so only the upper bound remains reachable.
#[test]
fn out_of_range_naca_code_is_rejected() {
    assert!(Airfoil::naca4(10_000, POINTS).is_err());
    assert!(Airfoil::naca4_str("-1", POINTS).is_err());
}

/// Ports `AirfoilExtraTest.NacaZeroThickness`.
#[test]
fn naca_zero_code_is_a_valid_zero_thickness_section() {
    let airfoil = Airfoil::naca4(0, POINTS).unwrap();
    assert!(airfoil.is_valid());
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.max_camber(), 0.0, epsilon = TOL);
}

/// Ports `AirfoilExtraTest.EvaluateThicknessNaca0012`.
#[test]
fn naca_0012_thickness_vanishes_at_the_leading_edge() {
    let airfoil = Airfoil::naca4(12, POINTS).unwrap();
    assert!(airfoil.thickness_at(0.3) > 0.0);
    assert_abs_diff_eq!(airfoil.thickness_at(0.0), 0.0, epsilon = 0.01);
}

/// Ports `AirfoilExtraTest.EvaluateCamberClamping`.
///
/// Queries outside `[0, 1]` clamp to the endpoints rather than extrapolating,
/// which keeps a mesher that overshoots by a rounding error from producing
/// nonsense camber.
#[test]
fn camber_evaluation_clamps_outside_the_chord() {
    let airfoil = Airfoil::naca4(2412, POINTS).unwrap();
    assert_abs_diff_eq!(airfoil.camber(-0.5), airfoil.camber(0.0), epsilon = TOL);
    assert_abs_diff_eq!(airfoil.camber(1.5), airfoil.camber(1.0), epsilon = TOL);
}

/// Ports `AirfoilExtraTest.InvalidAirfoilReturnsZero`.
#[test]
fn an_invalid_airfoil_evaluates_to_zero() {
    let airfoil = Airfoil::default();
    assert!(!airfoil.is_valid());
    assert_abs_diff_eq!(airfoil.camber(0.5), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.camber_slope(0.5), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(airfoil.thickness_at(0.5), 0.0, epsilon = TOL);
}

/// Ports `AirfoilExtraTest.SetName`.
#[test]
fn airfoil_name_can_be_overridden() {
    let mut airfoil = Airfoil::flat_plate();
    airfoil.set_name("Custom");
    assert_eq!(airfoil.name(), "Custom");
}

/// Ports `AirfoilExtraTest.CreateFromCoordinatesCamberLine`.
#[test]
fn airfoil_from_a_camber_line() {
    let airfoil = Airfoil::from_coordinates(&[0.0, 0.5, 1.0], &[0.0, 0.05, 0.0], true).unwrap();
    assert!(airfoil.is_valid());
    assert!(airfoil.max_camber() > 0.0);
    assert_abs_diff_eq!(airfoil.max_camber_position(), 0.5, epsilon = TOL);
}

/// Ports `AirfoilExtraTest.CreateFromCoordinatesMismatch` and
/// `CreateFromCoordinatesTooFew`.
#[test]
fn malformed_coordinate_input_is_rejected() {
    assert!(Airfoil::from_coordinates(&[0.0, 0.5, 1.0], &[0.0, 0.05], true).is_err());
    assert!(Airfoil::from_coordinates(&[0.5], &[0.0], true).is_err());
}

/// Ports `AirfoilExtraTest.CreateFromCoordinatesFullAirfoil`.
///
/// A closed loop traversed from the trailing edge round to the trailing edge
/// splits into camber and thickness. A symmetric loop has thickness but no
/// camber, which is what distinguishes this from the camber-line case.
#[test]
fn airfoil_from_a_full_closed_loop() {
    let airfoil = Airfoil::from_coordinates(
        &[1.0, 0.5, 0.0, 0.5, 1.0],
        &[0.0, 0.05, 0.0, -0.05, 0.0],
        false,
    )
    .unwrap();
    assert!(airfoil.is_valid());
    assert!(airfoil.max_thickness() > 0.0);
    assert_abs_diff_eq!(airfoil.max_camber(), 0.0, epsilon = 0.01);
}

/// Ports `AirfoilExtraTest.CreateFromCoordinateFileNotFound`.
#[test]
fn a_missing_coordinate_file_is_an_error() {
    assert!(Airfoil::from_coordinate_file("nonexistent.dat").is_err());
}

/// Ports `AirfoilExtraTest.CamberXY`.
#[test]
fn camber_arrays_span_the_unit_chord() {
    let airfoil = Airfoil::naca4(2412, POINTS).unwrap();
    let x = airfoil.camber_x();
    let y = airfoil.camber_y();
    assert_eq!(x.len(), y.len());
    assert!(x.len() > 2);
    assert_abs_diff_eq!(x[0], 0.0, epsilon = TOL);
    assert_abs_diff_eq!(x[x.len() - 1], 1.0, epsilon = TOL);
}

// ============================================================================
// Section
// ============================================================================

/// Ports `SectionTest.DefaultConstruction`.
#[test]
fn section_geometry_defaults() {
    let section = Section::default();
    assert_abs_diff_eq!(section.chord(), 1.0);
    assert_abs_diff_eq!(section.twist_degrees(), 0.0);
    assert_abs_diff_eq!(section.cl_alpha_factor(), 1.0, epsilon = TOL);
    assert!(!section.has_airfoil());
}

/// Ports `SectionTest.ParameterizedConstruction`.
#[test]
fn section_stores_degrees_and_reports_radians() {
    let position = Vec3::new(1.0, 2.0, 0.5);
    let section = Section::new(position, 0.3, 5.0);
    assert_abs_diff_eq!(*section.leading_edge(), position);
    assert_abs_diff_eq!(section.chord(), 0.3);
    assert_abs_diff_eq!(section.twist_degrees(), 5.0);
    assert_abs_diff_eq!(
        section.twist_radians(),
        5.0 * DEGREES_TO_RADIANS,
        epsilon = TOL
    );
}

/// Ports `SectionTest.AirfoilAssignment`.
#[test]
fn section_airfoil_assignment_keeps_the_name() {
    let mut section = Section::default();
    section.set_airfoil(Airfoil::naca4(2412, POINTS).unwrap());
    assert!(section.has_airfoil());
    assert_eq!(section.airfoil().name(), "NACA 2412");
}

/// Ports `SectionTest.ControlSurface`.
#[test]
fn section_carries_a_named_control_surface() {
    let mut section = Section::default();
    section.add_control_surface(aileron());
    assert_eq!(section.control_surfaces().len(), 1);
    assert_eq!(section.control_surfaces()[0].name, "aileron");
}

// ============================================================================
// Surface
// ============================================================================

/// Ports `SurfaceTest.BasicConstruction`.
#[test]
fn surface_basic_construction() {
    let surface = Surface::new("Main Wing");
    assert_eq!(surface.name(), "Main Wing");
    assert_eq!(surface.sections().len(), 0);
    assert!(surface.has_wake());
    assert!(surface.contributes_to_loads());
}

/// Ports `SurfaceTest.AddSections`.
#[test]
fn surface_sections_are_appended() {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 0.8, -2.0));
    assert_eq!(surface.sections().len(), 2);
}

/// Ports `SurfaceTest.PlanformAreaRectangular` and `PlanformAreaTapered`.
///
/// The tapered case is the one with teeth: trapezoidal integration gives
/// `(2 + 1)/2 * 4 = 6`, where a mean-chord shortcut would give something else.
#[test]
fn surface_planform_area_integrates_the_taper() {
    let mut rectangular = Surface::new("Wing");
    rectangular.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    rectangular.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(rectangular.planform_area(), 4.0, epsilon = TOL);

    let mut tapered = Surface::new("Wing");
    tapered.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 2.0, 0.0));
    tapered.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(tapered.planform_area(), 6.0, epsilon = TOL);
}

/// Ports `SurfaceTest.YDuplicateDoublesArea`.
#[test]
fn y_duplicate_doubles_the_planform_area() {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(surface.planform_area(), 2.0, epsilon = TOL);
    surface.set_y_duplicate(0.0);
    assert_abs_diff_eq!(surface.planform_area(), 4.0, epsilon = TOL);
}

/// Ports the four `MeanAerodynamicChord*` cases.
///
/// For a tapered wing the MAC is chord-weighted, so it sits strictly above the
/// arithmetic mean of 1.5 — the check that catches a plain average.
#[test]
fn surface_mean_aerodynamic_chord() {
    let mut rectangular = Surface::new("Wing");
    rectangular.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    rectangular.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(rectangular.mean_aerodynamic_chord(), 1.0, epsilon = TOL);

    let mut tapered = Surface::new("Wing");
    tapered.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 2.0, 0.0));
    tapered.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    let mac = tapered.mean_aerodynamic_chord();
    assert!(mac > 1.5);
    assert!(mac < 2.0);

    let mut single = Surface::new("Wing");
    single.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.5, 0.0));
    assert_abs_diff_eq!(single.mean_aerodynamic_chord(), 1.5, epsilon = TOL);

    assert_abs_diff_eq!(
        Surface::new("Wing").mean_aerodynamic_chord(),
        0.0,
        epsilon = TOL
    );
}

/// Ports `ComputeSpan`, `ComputeSpanYDuplicate` and `ComputeSpanEmpty`.
#[test]
fn surface_span_accounts_for_duplication() {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(surface.span(), 5.0, epsilon = TOL);
    surface.set_y_duplicate(0.0);
    assert_abs_diff_eq!(surface.span(), 10.0, epsilon = TOL);
    assert_abs_diff_eq!(Surface::new("Wing").span(), 0.0, epsilon = TOL);
}

/// Ports `PlanformAreaSingleSection`.
///
/// One section spans nothing, so it encloses no area — the degenerate case a
/// loop over intervals must handle by producing zero rather than indexing off
/// the end.
#[test]
fn a_single_section_encloses_no_area() {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    assert_abs_diff_eq!(surface.planform_area(), 0.0, epsilon = TOL);
}

/// Ports `HasViscousDragDefault` and `HasViscousDragNonZero`.
#[test]
fn surface_viscous_drag_flag_follows_the_polar() {
    let mut surface = Surface::new("Wing");
    assert!(!surface.has_viscous_drag());
    surface.set_cl_cd_polar([0.01, 0.0, 0.1, 0.0, 0.0, 0.0]);
    assert!(surface.has_viscous_drag());
    assert_abs_diff_eq!(surface.cl_cd_polar()[0], 0.01, epsilon = TOL);
}

/// Ports `ChordwiseSpacingNumeric` and `SpanwiseSpacingNumeric`.
#[test]
fn surface_numeric_spacing_follows_the_avl_sign_convention() {
    let mut surface = Surface::new("Wing");
    for (value, expected) in [
        (0.0, SpacingType::Uniform),
        (1.0, SpacingType::Cosine),
        (-1.0, SpacingType::Sine),
    ] {
        surface.set_numeric_chordwise_spacing(value);
        assert_eq!(surface.chordwise_spacing(), expected);
        surface.set_numeric_spanwise_spacing(value);
        assert_eq!(surface.spanwise_spacing(), expected);
    }
}

/// Ports `ClearYDuplicate`.
#[test]
fn y_duplicate_can_be_cleared() {
    let mut surface = Surface::new("Wing");
    surface.set_y_duplicate(0.0);
    assert!(surface.y_duplicate().is_some());
    surface.clear_y_duplicate();
    assert!(surface.y_duplicate().is_none());
}

/// Ports `AerodynamicOptions`.
#[test]
fn surface_aerodynamic_flags_round_trip() {
    let mut surface = Surface::new("Wing");
    surface.set_has_wake(false);
    assert!(!surface.has_wake());
    surface.set_contributes_to_loads(false);
    assert!(!surface.contributes_to_loads());
    surface.set_affected_by_alpha_beta(false);
    assert!(!surface.affected_by_alpha_beta());
}

/// Ports `CoordinateTransforms`.
#[test]
fn surface_transform_round_trips() {
    let mut surface = Surface::new("Wing");
    surface.set_angle_offset_degrees(5.0);
    assert_abs_diff_eq!(surface.angle_offset_degrees(), 5.0, epsilon = TOL);
    surface.set_translation(Vec3::new(1.0, 2.0, 3.0));
    assert_abs_diff_eq!(*surface.translation(), Vec3::new(1.0, 2.0, 3.0));
    surface.set_scale(Vec3::new(2.0, 2.0, 2.0));
    assert_abs_diff_eq!(*surface.scale(), Vec3::new(2.0, 2.0, 2.0));
}

/// Ports `ComponentIndex`.
#[test]
fn surface_component_index_round_trips() {
    let mut surface = Surface::new("Wing");
    surface.set_component_index(3);
    assert_eq!(surface.component_index(), 3);
}

/// Ports `SectionAccessors`.
#[test]
fn surface_section_accessors_read_and_mutate() {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 0.5, 0.0));
    assert_eq!(surface.sections().len(), 2);
    assert_abs_diff_eq!(surface.section(0).unwrap().chord(), 1.0, epsilon = TOL);
    assert_abs_diff_eq!(surface.section(1).unwrap().chord(), 0.5, epsilon = TOL);
    // Out-of-range access is an `Option`, not undefined behaviour.
    assert!(surface.section(2).is_none());

    surface.section_mut(0).unwrap().set_chord(0.8);
    assert_abs_diff_eq!(surface.section(0).unwrap().chord(), 0.8, epsilon = TOL);
    assert_eq!(surface.sections_mut().len(), 2);
}

// ============================================================================
// Aircraft
// ============================================================================

/// Ports `AircraftTest.DefaultValues`.
#[test]
fn aircraft_default_values() {
    let aircraft = Aircraft::default();
    assert_eq!(aircraft.surfaces().len(), 0);
    assert_eq!(aircraft.bodies().len(), 0);
    assert_abs_diff_eq!(aircraft.mach_number(), 0.0, epsilon = TOL);
}

/// Ports `AircraftTest.AddSurface` and `MoveSurface`.
#[test]
fn aircraft_accepts_surfaces_and_keeps_its_title() {
    let mut aircraft = Aircraft::new("Test Aircraft");
    let mut wing = Surface::new("Wing");
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    aircraft.add_surface(wing);
    assert_eq!(aircraft.surfaces().len(), 1);
    assert_eq!(aircraft.title(), "Test Aircraft");
    assert_eq!(aircraft.surfaces()[0].name(), "Wing");
}

/// Ports `AircraftTest.FindSurfaceByName`.
#[test]
fn aircraft_finds_a_surface_by_name() {
    let mut aircraft = Aircraft::default();
    let mut wing = Surface::new("Wing");
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    let mut tail = Surface::new("Tail");
    tail.add_section(Section::new(Vec3::new(3.0, 0.0, 0.0), 0.5, 0.0));
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    assert_eq!(aircraft.surface_named("Tail").unwrap().name(), "Tail");
    assert!(aircraft.surface_named("Missing").is_none());
}

/// Ports `AircraftTest.ReferenceQuantities`.
#[test]
fn aircraft_reference_quantities_round_trip() {
    let mut aircraft = Aircraft::default();
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.5);
    aircraft.set_reference_span(8.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    assert_abs_diff_eq!(aircraft.reference_area(), 10.0, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.reference_chord(), 1.5, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.reference_span(), 8.0, epsilon = TOL);
    assert_abs_diff_eq!(
        *aircraft.moment_reference_point(),
        Vec3::new(0.25, 0.0, 0.0)
    );
}

/// Ports `AircraftTest.FlowConditions`.
#[test]
fn aircraft_flow_conditions_round_trip() {
    let mut aircraft = Aircraft::default();
    aircraft.set_mach_number(0.3);
    assert_abs_diff_eq!(aircraft.mach_number(), 0.3, epsilon = TOL);
    aircraft.set_profile_drag_coefficient(0.01);
    assert_abs_diff_eq!(aircraft.profile_drag_coefficient(), 0.01, epsilon = TOL);
}

/// Ports `AircraftTest.Symmetry`.
#[test]
fn aircraft_symmetry_round_trips() {
    let mut aircraft = Aircraft::default();
    aircraft.set_y_symmetry(SymmetryType::Symmetric);
    assert_eq!(aircraft.y_symmetry(), SymmetryType::Symmetric);
    aircraft.set_z_symmetry(SymmetryType::Antisymmetric);
    assert_eq!(aircraft.z_symmetry(), SymmetryType::Antisymmetric);
    aircraft.set_z_symmetry_plane(1.0);
    assert_abs_diff_eq!(aircraft.z_symmetry_plane(), 1.0, epsilon = TOL);
}

/// Ports `AircraftTest.Bodies`.
#[test]
fn aircraft_accepts_bodies() {
    let mut aircraft = Aircraft::default();
    aircraft.add_body(Body::new("Fuselage"));
    assert_eq!(aircraft.bodies().len(), 1);
    assert_eq!(aircraft.bodies()[0].name(), "Fuselage");
}

/// Ports `AircraftTest.SurfaceAccessors`.
#[test]
fn aircraft_surface_accessors_read_and_mutate() {
    let mut aircraft = Aircraft::default();
    let mut wing = Surface::new("Wing");
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    aircraft.add_surface(wing);
    assert_eq!(aircraft.surfaces().len(), 1);
    assert_eq!(aircraft.surfaces()[0].name(), "Wing");
    aircraft.surfaces_mut()[0].set_name("Modified Wing");
    assert_eq!(aircraft.surfaces()[0].name(), "Modified Wing");
}

/// Ports `AircraftTest.ControlSurfaceNames`.
///
/// The same control appears on both sections; the aircraft-level list is a set
/// of names, so it must collapse to one entry rather than two.
#[test]
fn aircraft_control_surface_names_are_deduplicated() {
    let mut aircraft = Aircraft::default();
    let mut wing = Surface::new("Wing");
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0);
    root.add_control_surface(aileron());
    wing.add_section(root);
    let mut tip = Section::new(Vec3::new(0.0, 4.0, 0.0), 0.8, 0.0);
    tip.add_control_surface(aileron());
    wing.add_section(tip);
    aircraft.add_surface(wing);

    let names = aircraft.control_surface_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "aileron");
}

/// Ports `AircraftTest.DesignVariableNames`.
#[test]
fn aircraft_design_variable_names() {
    let mut aircraft = Aircraft::default();
    let mut wing = Surface::new("Wing");
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0);
    root.add_design_variable(DesignVariableDefinition {
        name: "twist".into(),
        weight: 1.0,
    });
    wing.add_section(root);
    aircraft.add_surface(wing);

    let names = aircraft.design_variable_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "twist");
}

/// Ports `ComputeTotalVortexCount`, `...YDuplicate`, `ComputeTotalStripCount`
/// and `...YDuplicate`.
///
/// Duplication doubles both counts. These figures have to agree with what the
/// mesher actually builds, not with a surface-level multiplication — that
/// disagreement was a real defect on this branch.
#[test]
fn aircraft_lattice_counts_include_duplication() {
    let build = |chordwise: usize, spanwise: usize, duplicated: bool| {
        let mut wing = Surface::new("Wing");
        wing.set_chordwise_panel_count(chordwise);
        wing.set_spanwise_panel_count(spanwise);
        if duplicated {
            wing.set_y_duplicate(0.0);
        }
        wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
        wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));
        let mut aircraft = Aircraft::default();
        aircraft.add_surface(wing);
        aircraft
    };
    assert_eq!(build(4, 6, false).total_vortex_count(), 24);
    assert_eq!(build(4, 6, true).total_vortex_count(), 48);
    assert_eq!(build(1, 8, false).total_strip_count(), 8);
    assert_eq!(build(1, 8, true).total_strip_count(), 16);
}

/// Ports `ComputeTotalPlanformArea` and `ComputeTotalSpan`.
///
/// Area sums across surfaces (5.0 + 0.6); span is the maximum, not the sum,
/// since surfaces overlap spanwise.
#[test]
fn aircraft_totals_sum_area_but_take_the_maximum_span() {
    let mut wing = Surface::new("Wing");
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));
    let mut tail = Surface::new("Tail");
    tail.add_section(Section::new(Vec3::new(4.0, 0.0, 0.0), 0.5, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.3, 0.0));

    let mut aircraft = Aircraft::default();
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    assert_abs_diff_eq!(aircraft.total_planform_area(), 5.6, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.total_span(), 5.0, epsilon = TOL);
}

// ============================================================================
// AVL reader
// ============================================================================

/// Ports `AvlReaderTest.ParseSimpleGeometry`.
#[test]
fn avl_parses_simple_geometry() {
    let content = "\
Test Wing
0.0
0 0 0
10.0 1.0 10.0
0.5 0.0 0.0
0.0

SURFACE
Wing
12 1.0 20 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.title(), "Test Wing");
    assert_abs_diff_eq!(aircraft.mach_number(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.reference_area(), 10.0, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.reference_chord(), 1.0, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.reference_span(), 10.0, epsilon = TOL);
    assert_eq!(aircraft.surfaces().len(), 1);
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.name(), "Wing");
    assert_eq!(wing.chordwise_panel_count(), 12);
    assert_eq!(wing.spanwise_panel_count(), 20);
    assert_eq!(wing.sections().len(), 2);
}

/// Ports `AvlReaderTest.ParseWithNacaAirfoil`.
#[test]
fn avl_parses_a_naca_keyword() {
    let content = "\
NACA Test
0.1
0 0 0
1.0 0.25 4.0
0.0 0.0 0.0
0.0

SURFACE
Wing
8 1.0

SECTION
0.0 0.0 0.0 0.5 0.0
NACA
2412

SECTION
0.0 2.0 0.0 0.5 0.0
NACA
2412
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.sections().len(), 2);
    assert!(wing.sections()[0].has_airfoil());
    assert_eq!(wing.sections()[0].airfoil().name(), "NACA 2412");
}

/// Ports `AvlReaderTest.ParseWithYDuplicate`.
#[test]
fn avl_parses_y_duplicate() {
    let content = "\
Symmetric Wing
0.0
0 0 0
2.0 0.5 4.0
0.0 0.0 0.0

SURFACE
Wing
6 1.0
YDUPLICATE
0.0

SECTION
0.0 0.0 0.0 0.5 0.0

SECTION
0.0 2.0 0.0 0.5 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let wing = &aircraft.surfaces()[0];
    assert_abs_diff_eq!(wing.y_duplicate().unwrap(), 0.0, epsilon = TOL);
}

/// Ports `AvlReaderTest.ParseWithControlSurface`.
#[test]
fn avl_parses_a_control_surface_on_both_sections() {
    let content = "\
Control Test
0.0
0 0 0
1.0 0.25 2.0
0.0 0.0 0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 0.5 0.0
CONTROL
aileron 1.0 0.75 0 1 0 1.0

SECTION
0.0 1.0 0.0 0.5 0.0
CONTROL
aileron 1.0 0.75 0 1 0 -1.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let names = aircraft.control_surface_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "aileron");
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.sections()[0].control_surfaces().len(), 1);
    assert_eq!(wing.sections()[0].control_surfaces()[0].name, "aileron");
}

/// Ports `AvlReaderExtraTest.ReadFileNotFound` and `ParseEmptyContent`.
#[test]
fn avl_rejects_a_missing_file_and_empty_content() {
    assert!(AvlReader::default().read_path("nonexistent.avl").is_err());
    assert!(AvlReader::default().read_str("").is_err());
}

/// Ports `AvlReaderExtraTest.ParseWithComments`.
///
/// Both comment markers have to be stripped before the title is read, or the
/// first comment becomes the aircraft's name.
#[test]
fn avl_strips_both_comment_markers() {
    let content = "\
# This is a comment
! Another comment
Test Wing
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.title(), "Test Wing");
}

/// Ports `AvlReaderExtraTest.ParseWithProfileDrag`.
#[test]
fn avl_parses_the_profile_drag_coefficient() {
    let content = "\
Wing with drag
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.02

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_abs_diff_eq!(aircraft.profile_drag_coefficient(), 0.02, epsilon = TOL);
}

/// Ports `AvlReaderExtraTest.ParseWithSymmetry`.
#[test]
fn avl_parses_the_symmetry_line() {
    let content = "\
Symmetric Wing
0.0
1 -1 0.5
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.y_symmetry(), SymmetryType::Symmetric);
    assert_eq!(aircraft.z_symmetry(), SymmetryType::Antisymmetric);
    assert_abs_diff_eq!(aircraft.z_symmetry_plane(), 0.5, epsilon = TOL);
}

/// Ports `AvlReaderExtraTest.ParseSurfaceModifiers`.
///
/// Every surface-level keyword in one file, including the three negating flags
/// and the `CDCL` polar.
#[test]
fn avl_parses_every_surface_modifier() {
    let content = "\
Modifiers Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0 8 1.0
SCALE
2.0 2.0 2.0
TRANSLATE
1.0 0.0 0.0
ANGLE
5.0
COMPONENT
2
NOWAKE
NOLOAD
NOALBE
CDCL
0.01 0.0 0.005 0.0 0.0 0.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let wing = &aircraft.surfaces()[0];
    assert_abs_diff_eq!(*wing.scale(), Vec3::new(2.0, 2.0, 2.0));
    assert_abs_diff_eq!(*wing.translation(), Vec3::new(1.0, 0.0, 0.0));
    assert_abs_diff_eq!(wing.angle_offset_degrees(), 5.0, epsilon = TOL);
    assert_eq!(wing.component_index(), 2);
    assert!(!wing.has_wake());
    assert!(!wing.contributes_to_loads());
    assert!(!wing.affected_by_alpha_beta());
    assert!(wing.has_viscous_drag());
    assert_eq!(wing.spanwise_panel_count(), 8);
}

/// Ports `AvlReaderExtraTest.ParseSectionSpanwise`.
#[test]
fn avl_parses_section_level_spanwise_and_claf() {
    let content = "\
Section Spanwise
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 3.0 4 -1.0
CLAF
1.5

SECTION
0.0 5.0 0.0 0.8 -2.0 4 1.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.sections().len(), 2);
    assert_abs_diff_eq!(wing.sections()[0].twist_degrees(), 3.0, epsilon = TOL);
    assert_abs_diff_eq!(wing.sections()[0].cl_alpha_factor(), 1.5, epsilon = TOL);
    assert_abs_diff_eq!(wing.sections()[1].chord(), 0.8, epsilon = TOL);
}

/// Ports `AvlReaderExtraTest.ParseWithDesignVariable`.
#[test]
fn avl_parses_design_variables() {
    let content = "\
Design Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0
DESIGN
twist 1.0

SECTION
0.0 5.0 0.0 1.0 0.0
DESIGN
twist 0.5
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    let names = aircraft.design_variable_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "twist");
}

/// Ports `AvlReaderExtraTest.ParseWithBody`.
#[test]
fn avl_parses_a_body_with_modifiers() {
    let content = "\
Body Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

BODY
Fuselage
20 1.0
SCALE
1.0 1.0 1.0
TRANSLATE
0.0 0.0 0.0
YDUPLICATE
BSEC
0.0 0.0 0.0 0.2 0.2
BSEC
2.0 0.0 0.0 0.5 0.5
BSEC
5.0 0.0 0.0 0.1 0.1

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.bodies().len(), 1);
    assert_eq!(aircraft.bodies()[0].name(), "Fuselage");
    assert!(aircraft.bodies()[0].is_y_duplicated());
    assert_eq!(aircraft.surfaces().len(), 1);
}

/// Ports `AvlReaderExtraTest.ParseBodyAtEnd`.
///
/// A body that runs to end-of-file has no following keyword to terminate it,
/// so it has to be flushed when the input is exhausted.
#[test]
fn avl_flushes_a_body_at_end_of_file() {
    let content = "\
Body at End
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

BODY
Fuselage
20 1.0
BSEC
0.0 0.0 0.0 0.2 0.2
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.bodies().len(), 1);
    assert_eq!(aircraft.bodies()[0].name(), "Fuselage");
}

/// Ports `AvlReaderExtraTest.ParseVerbose`.
#[test]
fn avl_verbose_mode_parses_identically() {
    let content = "\
Verbose Test
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let mut reader = AvlReader::default();
    reader.set_verbose(true);
    let aircraft = reader.read_str(content).unwrap();
    assert_eq!(aircraft.title(), "Verbose Test");
}

/// Ports `AvlReaderExtraTest.ParseMomentReference`.
#[test]
fn avl_parses_mach_and_the_moment_reference() {
    let content = "\
MomRef Test
0.5
0 0 0
10.0 1.0 10.0
0.25 0.0 0.1
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_abs_diff_eq!(aircraft.mach_number(), 0.5, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.moment_reference_point().x, 0.25, epsilon = TOL);
    assert_abs_diff_eq!(aircraft.moment_reference_point().z, 0.1, epsilon = TOL);
}

/// Ports `AvlReaderExtraTest.ParseMultipleSurfaces`.
#[test]
fn avl_parses_multiple_surfaces_in_order() {
    let content = "\
Multi Surface
0.0
0 0 0
10.0 1.0 10.0
0.0 0.0 0.0
0.0

SURFACE
Wing
6 1.0

SECTION
0.0 0.0 0.0 1.0 0.0

SECTION
0.0 5.0 0.0 1.0 0.0

SURFACE
Tail
4 1.0

SECTION
4.0 0.0 0.0 0.5 0.0

SECTION
4.0 1.5 0.0 0.3 0.0
";
    let aircraft = AvlReader::default().read_str(content).unwrap();
    assert_eq!(aircraft.surfaces().len(), 2);
    assert_eq!(aircraft.surfaces()[0].name(), "Wing");
    assert_eq!(aircraft.surfaces()[1].name(), "Tail");
}
