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

//! Ported from the frozen v0.2.0 `section_test.cpp`.
//!
//! `Section` is the per-station geometry record. Two things here are easy to
//! get wrong and are asserted directly: twist is stored in degrees but
//! consumed in radians, and the numeric spanwise-spacing setter follows AVL's
//! sign convention rather than a magnitude.

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::{PI, SpacingType, Vec3};
use zephyr_vlm::geometry::{Airfoil, ControlSurfaceDefinition, DesignVariableDefinition, Section};

const TOL: f64 = 1.0e-10;

/// Ports `DefaultConstructor`.
#[test]
fn section_default_constructor() {
    let section = Section::default();
    assert_abs_diff_eq!(section.chord(), 1.0, epsilon = TOL);
    assert_abs_diff_eq!(section.twist_degrees(), 0.0, epsilon = TOL);
    assert!(!section.has_airfoil());
}

/// Ports `ParameterizedConstructor`.
#[test]
fn section_parameterized_constructor() {
    let section = Section::new(Vec3::new(1.0, 2.0, 3.0), 0.5, 5.0);
    assert_abs_diff_eq!(*section.leading_edge(), Vec3::new(1.0, 2.0, 3.0));
    assert_abs_diff_eq!(section.chord(), 0.5, epsilon = TOL);
    assert_abs_diff_eq!(section.twist_degrees(), 5.0, epsilon = TOL);
    assert!(!section.has_airfoil());
}

/// Ports `ConstructorWithAirfoil` and `GetAirfoilWithCustomReturnsCustom`.
#[test]
fn section_constructed_with_an_airfoil_keeps_it() {
    let section = Section::with_airfoil(
        Vec3::zeros(),
        1.0,
        0.0,
        Airfoil::naca4_str("2412", 61).unwrap(),
    );
    assert!(section.has_airfoil());
    assert!(section.airfoil().thickness_at(0.5) > 0.0);
}

/// Ports `GetAirfoilWithoutCustomReturnsFlat`.
///
/// Reading the airfoil of a section that has none must hand back a flat plate
/// rather than fail, so the mesher can treat every section uniformly.
#[test]
fn section_without_an_airfoil_reads_as_a_flat_plate() {
    let section = Section::new(Vec3::zeros(), 1.0, 0.0);
    assert_abs_diff_eq!(section.airfoil().thickness_at(0.5), 0.0, epsilon = TOL);
}

/// Ports `GetAirfoilMutableCreatesFlat`.
///
/// The mutable accessor materializes the implicit flat plate, so the section
/// reports as having an airfoil afterwards even though nothing was assigned.
#[test]
fn section_mutable_airfoil_accessor_materializes_a_flat_plate() {
    let mut section = Section::new(Vec3::zeros(), 1.0, 0.0);
    assert!(!section.has_airfoil());
    let thickness = section.airfoil_mut().thickness_at(0.5);
    assert!(section.has_airfoil());
    assert_abs_diff_eq!(thickness, 0.0, epsilon = TOL);
}

/// Ports `SetAirfoil`.
#[test]
fn section_set_airfoil() {
    let mut section = Section::default();
    assert!(!section.has_airfoil());
    section.set_airfoil(Airfoil::naca4_str("0012", 61).unwrap());
    assert!(section.has_airfoil());
}

/// Ports `SetLeadingEdgePosition`, `SetChord` and `SetTwistAngle`.
#[test]
fn section_geometry_setters_round_trip() {
    let mut section = Section::default();
    section.set_leading_edge(Vec3::new(5.0, 10.0, -1.0));
    assert_abs_diff_eq!(*section.leading_edge(), Vec3::new(5.0, 10.0, -1.0));
    section.set_chord(2.5);
    assert_abs_diff_eq!(section.chord(), 2.5, epsilon = TOL);
    section.set_twist_degrees(3.0);
    assert_abs_diff_eq!(section.twist_degrees(), 3.0, epsilon = TOL);
}

/// Ports `TwistAngleRadians`.
///
/// The stored value is degrees; the consumed value is radians. This is the
/// deliberate unit split, so a 90 degree twist must read back as pi/2.
#[test]
fn section_twist_converts_degrees_to_radians() {
    let mut section = Section::default();
    section.set_twist_degrees(90.0);
    assert_abs_diff_eq!(section.twist_radians(), PI / 2.0, epsilon = 1.0e-6);
}

/// Ports `ClAlphaFactor`.
#[test]
fn section_cl_alpha_factor_defaults_to_one() {
    let mut section = Section::default();
    assert_abs_diff_eq!(section.cl_alpha_factor(), 1.0, epsilon = TOL);
    section.set_cl_alpha_factor(0.8);
    assert_abs_diff_eq!(section.cl_alpha_factor(), 0.8, epsilon = TOL);
}

/// Ports `AddControlSurface`.
#[test]
fn section_control_surfaces_are_appended() {
    let mut section = Section::default();
    assert_eq!(section.control_surfaces().len(), 0);
    section.add_control_surface(ControlSurfaceDefinition {
        name: "aileron".into(),
        hinge_x_fraction: 0.75,
        hinge_vector: Vec3::new(0.0, 1.0, 0.0),
        sign: 1.0,
        gain: 1.0,
    });
    assert_eq!(section.control_surfaces().len(), 1);
    assert_eq!(section.control_surfaces()[0].name, "aileron");
    assert_abs_diff_eq!(
        section.control_surfaces()[0].hinge_x_fraction,
        0.75,
        epsilon = TOL
    );
}

/// Ports `AddDesignVariable`.
#[test]
fn section_design_variables_are_appended() {
    let mut section = Section::default();
    assert!(section.design_variables().is_empty());
    section.add_design_variable(DesignVariableDefinition {
        name: "twist".into(),
        weight: 1.5,
    });
    assert_eq!(section.design_variables().len(), 1);
    assert_eq!(section.design_variables()[0].name, "twist");
    assert_abs_diff_eq!(section.design_variables()[0].weight, 1.5, epsilon = TOL);
}

/// Ports `SpanwisePanelCount`.
///
/// Zero means "no override": the surface-level count applies instead.
#[test]
fn section_spanwise_panel_count_defaults_to_no_override() {
    let mut section = Section::default();
    assert_eq!(section.spanwise_panel_count(), 0);
    section.set_spanwise_panel_count(8);
    assert_eq!(section.spanwise_panel_count(), 8);
}

/// Ports `SpanwiseSpacingEnum`.
#[test]
fn section_spanwise_spacing_enum_round_trips() {
    let mut section = Section::default();
    assert_eq!(section.spanwise_spacing(), SpacingType::None);
    section.set_spanwise_spacing(SpacingType::Cosine);
    assert_eq!(section.spanwise_spacing(), SpacingType::Cosine);
}

/// Ports `SpanwiseSpacingNumericUniform`, `...Cosine` and `...Sine`.
///
/// AVL encodes the distribution in the *sign* of the spacing field, not its
/// magnitude: zero is uniform, positive is cosine, negative is sine.
#[test]
fn section_numeric_spanwise_spacing_follows_the_avl_sign_convention() {
    let mut section = Section::default();
    section.set_numeric_spanwise_spacing(0.0);
    assert_eq!(section.spanwise_spacing(), SpacingType::Uniform);
    section.set_numeric_spanwise_spacing(1.0);
    assert_eq!(section.spanwise_spacing(), SpacingType::Cosine);
    section.set_numeric_spanwise_spacing(-1.0);
    assert_eq!(section.spanwise_spacing(), SpacingType::Sine);
}

/// Ports `ViscousPolar`.
#[test]
fn section_viscous_polar_defaults_to_zero_and_round_trips() {
    let mut section = Section::default();
    for value in section.viscous_polar() {
        assert_abs_diff_eq!(*value, 0.0, epsilon = TOL);
    }
    section.set_viscous_polar([0.01, 0.0, 0.005, 0.0, 0.0, 0.0]);
    assert_abs_diff_eq!(section.viscous_polar()[0], 0.01, epsilon = TOL);
    assert_abs_diff_eq!(section.viscous_polar()[2], 0.005, epsilon = TOL);
}
