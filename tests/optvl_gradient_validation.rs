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

//! Ported from the frozen v0.2.0 `optvl_gradient_validation_test.cpp`.
//!
//! End-to-end validation of the adjoint against full re-solves on a
//! wing-plus-tail configuration with sweep, dihedral, taper and twist — the
//! most complete gradient check in the frozen suite.
//!
//! # A unit bug reproduced from the fixture
//!
//! The C++ built its wing tip with
//! `Section(Vec3(0.2, 4.0, 0.3), 0.8, -2.0 * kDegreesToRadians)`. `Section`
//! stores twist in *degrees*, so that argument is `-0.0349` degrees, not the
//! `-2` degrees it reads as. The fixture is therefore very nearly untwisted.
//!
//! It is reproduced exactly rather than corrected: every tolerance in this file
//! was chosen against that geometry, and quietly changing the twist by a factor
//! of 57 would alter what the assertions mean. The Rust spelling makes the
//! conversion explicit so the mistake is visible instead of hidden in a
//! constant.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{AeroOutput, GeometryVariable, GradientCalculator};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

const FD_STEP: f64 = 1.0e-6;
/// The frozen 5% relative band and 1e-4 absolute floor.
const RELATIVE_TOLERANCE: f64 = 0.05;
const ABSOLUTE_TOLERANCE: f64 = 1.0e-4;

/// The frozen fixture: a swept, dihedralled, tapered wing plus a horizontal
/// tail aft of the moment reference.
fn aircraft() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(8);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    // See the module docs: the frozen fixture passes a radian value into a
    // degree field, giving a near-zero tip twist. Reproduced verbatim.
    wing.add_section(Section::new(
        Vec3::new(0.2, 4.0, 0.3),
        0.8,
        -2.0 * DEGREES_TO_RADIANS,
    ));

    let mut htail = Surface::new("HTail");
    htail.set_chordwise_panel_count(3);
    htail.set_spanwise_panel_count(4);
    htail.add_section(Section::new(Vec3::new(4.0, 0.0, 0.2), 0.5, 0.0));
    htail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.25), 0.4, 0.0));

    let mut aircraft = Aircraft::new("ValidationAircraft");
    aircraft.add_surface(wing);
    aircraft.add_surface(htail);
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(8.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft
}

fn condition(alpha_degrees: f64, mach: f64) -> FlightCondition {
    FlightCondition {
        alpha: alpha_degrees * DEGREES_TO_RADIANS,
        mach,
        ..FlightCondition::default()
    }
}

fn calculator() -> GradientCalculator {
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft());
    calculator
}

/// One coefficient from a full re-solve of the given aircraft and condition.
fn coefficient(aircraft: &Aircraft, condition: &FlightCondition, output: AeroOutput) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let solution = solver.solve(condition).unwrap();
    let results =
        ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, condition);
    match output {
        AeroOutput::Lift => results.coefficients.lift,
        AeroOutput::Drag => results.coefficients.drag,
        AeroOutput::Side => results.coefficients.side,
        AeroOutput::Roll => results.coefficients.roll,
        AeroOutput::Pitch => results.coefficients.pitch,
        AeroOutput::Yaw => results.coefficients.yaw,
        AeroOutput::InducedDrag => results.far_field_induced_drag,
    }
}

/// Centered difference of one coefficient with respect to one flight variable.
fn difference(
    base: &FlightCondition,
    output: AeroOutput,
    step: f64,
    perturb: impl Fn(&mut FlightCondition, f64),
) -> f64 {
    let aircraft = aircraft();
    let mut plus = base.clone();
    let mut minus = base.clone();
    perturb(&mut plus, step);
    perturb(&mut minus, -step);
    (coefficient(&aircraft, &plus, output) - coefficient(&aircraft, &minus, output)) / (2.0 * step)
}

/// The frozen comparison: 5% relative with a 1e-4 absolute floor.
fn assert_matches(adjoint: f64, finite: f64) {
    assert_abs_diff_eq!(
        adjoint,
        finite,
        epsilon = finite.abs() * RELATIVE_TOLERANCE + ABSOLUTE_TOLERANCE
    );
}

/// Ports `AlphaDerivativesMatchFD`.
#[test]
fn alpha_derivatives_match_finite_differences() {
    let base = condition(5.0, 0.5);
    let gradients = calculator().compute_alpha_gradients(&base).unwrap();

    let lift = difference(&base, AeroOutput::Lift, FD_STEP, |c, h| c.alpha += h);
    assert_matches(gradients.coefficients.lift, lift);

    let pitch = difference(&base, AeroOutput::Pitch, FD_STEP, |c, h| c.alpha += h);
    assert_matches(gradients.coefficients.pitch, pitch);

    // At positive incidence, both drag and induced drag rise with alpha.
    let drag = difference(&base, AeroOutput::Drag, FD_STEP, |c, h| c.alpha += h);
    assert!(gradients.coefficients.drag > 0.0);
    assert!(drag > 0.0);
    let induced = difference(&base, AeroOutput::InducedDrag, FD_STEP, |c, h| c.alpha += h);
    assert!(induced > 0.0);
}

/// Ports `BetaDerivativesMatchFD`.
#[test]
fn beta_derivatives_match_finite_differences() {
    let base = condition(5.0, 0.5);
    let gradients = calculator().compute_beta_gradients(&base).unwrap();

    let side = difference(&base, AeroOutput::Side, FD_STEP, |c, h| c.beta += h);
    assert_matches(gradients.coefficients.side, side);

    let roll = difference(&base, AeroOutput::Roll, FD_STEP, |c, h| c.beta += h);
    assert_matches(gradients.coefficients.roll, roll);

    let yaw = difference(&base, AeroOutput::Yaw, FD_STEP, |c, h| c.beta += h);
    assert!(gradients.coefficients.yaw.is_finite());
    assert!(yaw.is_finite());
}

/// Ports `MachDerivativesMatchFD`.
///
/// The frozen case uses a larger Mach step than the angle steps, because the
/// Prandtl-Glauert factor is much less sensitive per unit input.
#[test]
fn mach_derivatives_match_finite_differences() {
    let base = condition(5.0, 0.5);
    let gradients = calculator()
        .compute_flight_condition_gradients(&base)
        .unwrap();
    let step = 0.001;

    let lift = difference(&base, AeroOutput::Lift, step, |c, h| c.mach += h);
    assert_matches(gradients.by_mach.coefficients.lift, lift);

    let pitch = difference(&base, AeroOutput::Pitch, step, |c, h| c.mach += h);
    assert_matches(gradients.by_mach.coefficients.pitch, pitch);

    let induced = difference(&base, AeroOutput::InducedDrag, step, |c, h| c.mach += h);
    assert_matches(gradients.by_mach.induced_drag, induced);
}

/// Ports `FullGradientsConsistency`.
///
/// The batched flight-condition pass and the single-input methods must agree to
/// 1%, since they are two routes through the same adjoint.
#[test]
fn the_batched_pass_agrees_with_the_single_input_methods() {
    let base = condition(5.0, 0.3);
    let calculator = calculator();
    let full = calculator
        .compute_flight_condition_gradients(&base)
        .unwrap();
    let alpha = calculator.compute_alpha_gradients(&base).unwrap();
    let beta = calculator.compute_beta_gradients(&base).unwrap();

    assert_abs_diff_eq!(
        full.by_alpha.coefficients.lift,
        alpha.coefficients.lift,
        epsilon = alpha.coefficients.lift.abs() * 0.01
    );
    assert_abs_diff_eq!(
        full.by_beta.coefficients.side,
        beta.coefficients.side,
        epsilon = beta.coefficients.side.abs() * 0.01 + ABSOLUTE_TOLERANCE
    );
}

/// Ports `GeometryGradientsVsFiniteDifferences`.
#[test]
fn every_geometry_derivative_is_finite_on_the_validation_aircraft() {
    let gradients = calculator()
        .compute_geometry_gradients(&condition(5.0, 0.3))
        .unwrap();
    assert!(!gradients.sections.is_empty());
    for section in &gradients.sections {
        for group in [
            &section.by_chord,
            &section.by_twist,
            &section.by_leading_edge_x,
            &section.by_leading_edge_y,
            &section.by_leading_edge_z,
        ] {
            assert!(group.coefficients.lift.is_finite());
        }
    }
}

/// Perturbs one section design variable and returns the re-solved coefficient.
fn perturbed_coefficient(
    surface: usize,
    section: usize,
    variable: GeometryVariable,
    delta: f64,
    condition: &FlightCondition,
    output: AeroOutput,
) -> f64 {
    let mut perturbed = aircraft();
    {
        let target = perturbed.surfaces_mut()[surface]
            .section_mut(section)
            .unwrap();
        match variable {
            GeometryVariable::Chord => target.set_chord(target.chord() + delta),
            GeometryVariable::Twist => {
                target.set_twist_degrees(target.twist_degrees() + delta);
            }
            GeometryVariable::LeadingEdgeX => {
                let mut position = *target.leading_edge();
                position.x += delta;
                target.set_leading_edge(position);
            }
            GeometryVariable::LeadingEdgeY => {
                let mut position = *target.leading_edge();
                position.y += delta;
                target.set_leading_edge(position);
            }
            GeometryVariable::LeadingEdgeZ => {
                let mut position = *target.leading_edge();
                position.z += delta;
                target.set_leading_edge(position);
            }
        }
    }
    coefficient(&perturbed, condition, output)
}

/// Ports `AdjointGeometryGradientMatchesFD`.
///
/// The strongest case in the frozen suite: every section, every one of the five
/// design variables, against a full re-solve of the perturbed aircraft. Lift,
/// drag and pitching moment are all checked, so an error confined to one output
/// cannot hide.
#[test]
fn adjoint_geometry_gradients_match_resolves_for_every_design_variable() {
    let base = condition(5.0, 0.3);
    let gradients = calculator().compute_geometry_gradients(&base).unwrap();
    let aircraft = aircraft();

    // Section gradients are reported surface-major, matching the aircraft's own
    // surface and section order.
    let mut flat = Vec::new();
    for (surface_index, surface) in aircraft.surfaces().iter().enumerate() {
        for section_index in 0..surface.sections().len() {
            flat.push((surface_index, section_index));
        }
    }
    assert_eq!(flat.len(), gradients.sections.len());

    for (entry, (surface_index, section_index)) in gradients.sections.iter().zip(flat) {
        for (variable, group) in [
            (GeometryVariable::Chord, &entry.by_chord),
            (GeometryVariable::Twist, &entry.by_twist),
            (GeometryVariable::LeadingEdgeX, &entry.by_leading_edge_x),
            (GeometryVariable::LeadingEdgeY, &entry.by_leading_edge_y),
            (GeometryVariable::LeadingEdgeZ, &entry.by_leading_edge_z),
        ] {
            for (output, exact) in [
                (AeroOutput::Lift, group.coefficients.lift),
                (AeroOutput::Drag, group.coefficients.drag),
                (AeroOutput::Pitch, group.coefficients.pitch),
            ] {
                let step = 1.0e-6;
                let plus = perturbed_coefficient(
                    surface_index,
                    section_index,
                    variable,
                    step,
                    &base,
                    output,
                );
                let minus = perturbed_coefficient(
                    surface_index,
                    section_index,
                    variable,
                    -step,
                    &base,
                    output,
                );
                let finite = (plus - minus) / (2.0 * step);
                assert_abs_diff_eq!(
                    exact,
                    finite,
                    epsilon = finite.abs() * RELATIVE_TOLERANCE + ABSOLUTE_TOLERANCE
                );
            }
        }
    }
}

/// Ports `TwistGradientSignIsCorrect`.
#[test]
fn twisting_either_wing_station_increases_lift() {
    let gradients = calculator()
        .compute_geometry_gradients(&condition(5.0, 0.0))
        .unwrap();
    // The first two entries are the wing's root and tip.
    for section in &gradients.sections[..2] {
        assert!(
            section.by_twist.coefficients.lift > 0.0,
            "section {} gave {}",
            section.section_index,
            section.by_twist.coefficients.lift
        );
    }
}

/// Ports `TailTwistGradientAffectsCmStrongly`.
///
/// The tail sits aft of the moment reference, so twisting it pitches the
/// aircraft nose down — and by a margin well clear of numerical noise.
#[test]
fn tail_twist_drives_the_pitching_moment_strongly() {
    let gradients = calculator()
        .compute_geometry_gradients(&condition(5.0, 0.0))
        .unwrap();
    for section in &gradients.sections[2..] {
        assert_eq!(section.surface_name, "HTail");
        assert!(section.by_twist.coefficients.pitch < 0.0);
        assert!(section.by_twist.coefficients.pitch.abs() > 0.01);
    }
}

/// Ports `LiftCurveSlopeIsReasonable`.
#[test]
fn the_lift_curve_slope_matches_optvl() {
    let gradients = calculator()
        .compute_alpha_gradients(&condition(5.0, 0.3))
        .unwrap();
    let reference: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/avl_reference/optvl.json")).unwrap();
    let case = reference["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|case| case["name"] == "wing_tail")
        .unwrap();
    assert_abs_diff_eq!(
        gradients.coefficients.lift,
        case["alpha_derivatives"]["CL"].as_f64().unwrap(),
        epsilon = 2.0e-7
    );
}

/// Ports `RollDampingIsNegative` and `PitchRateDerivativesAreComputed`.
#[test]
fn rate_derivatives_have_stable_signs_and_are_populated() {
    let gradients = calculator()
        .compute_flight_condition_gradients(&condition(5.0, 0.0))
        .unwrap();
    assert!(gradients.by_roll_rate.coefficients.roll < 0.0);
    assert!(gradients.by_pitch_rate.coefficients.pitch != 0.0);
    assert!(gradients.by_pitch_rate.coefficients.lift != 0.0);
}

/// Ports `CompressibilityIncreasesCoefficients`.
#[test]
fn compressibility_steepens_the_lift_curve_slope() {
    let calculator = calculator();
    let low = calculator
        .compute_alpha_gradients(&condition(5.0, 0.1))
        .unwrap();
    let high = calculator
        .compute_alpha_gradients(&condition(5.0, 0.6))
        .unwrap();
    assert!(high.coefficients.lift > low.coefficients.lift);
}

/// Ports the symmetric-configuration case.
///
/// At zero sideslip a laterally symmetric aircraft produces almost no side
/// force or yawing moment in response to alpha; the frozen bounds are relative
/// to the corresponding longitudinal derivative rather than absolute.
#[test]
fn a_symmetric_aircraft_has_negligible_lateral_response_to_alpha() {
    let gradients = calculator()
        .compute_alpha_gradients(&condition(5.0, 0.0))
        .unwrap();
    assert!(gradients.coefficients.side.abs() < gradients.coefficients.lift.abs() * 0.1);
    assert!(gradients.coefficients.yaw.abs() < gradients.coefficients.pitch.abs() * 0.5);
}

/// Ports `ZeroAlphaGivesMinimumDragSlope`.
///
/// Induced drag is quadratic in lift, so its slope passes through zero where
/// the lift does.
#[test]
fn the_drag_slope_is_near_zero_at_zero_incidence() {
    let gradients = calculator()
        .compute_alpha_gradients(&condition(0.0, 0.3))
        .unwrap();
    assert_abs_diff_eq!(gradients.coefficients.drag, 0.0, epsilon = 0.05);
}
