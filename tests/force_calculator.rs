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

//! Ported from the frozen v0.2.0 `force_calculator_test.cpp` and
//! `stability_derivatives_test.cpp`.
//!
//! These are the physics checks rather than type contracts: sign conventions,
//! the quadratic lift-drag relation, agreement between the near-field sum and
//! the Trefftz-plane integral, and the derivative signs that make an aircraft
//! flyable.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::{AerodynamicResults, ForceCalculator, StabilityDerivativeCalculator};
use zephyr_vlm::core::{DEGREES_TO_RADIANS, PI, SpacingType, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

const TOL: f64 = 1.0e-6;

/// The frozen fixture: a flat rectangular wing centred on `y = 0`, cosine
/// spaced in both directions.
fn rectangular_wing(span: f64, chord: f64, chordwise: usize, spanwise: usize) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    wing.set_chordwise_spacing(SpacingType::Cosine);
    wing.set_spanwise_spacing(SpacingType::Cosine);
    wing.add_section(Section::new(Vec3::new(0.0, -span / 2.0, 0.0), chord, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));

    let mut aircraft = Aircraft::new("Test Wing");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn solve(aircraft: &Aircraft, alpha_degrees: f64) -> (VortexLatticeSolver, AerodynamicResults) {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let condition = FlightCondition {
        alpha: alpha_degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    assert!(solution.is_valid);
    let results =
        ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, &condition);
    (solver, results)
}

// ============================================================================
// force_calculator_test.cpp
// ============================================================================

/// Ports `ZeroAlphaGivesZeroForces`.
#[test]
fn zero_alpha_on_a_flat_wing_gives_zero_lift() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let (_, results) = solve(&aircraft, 0.0);
    assert!(results.is_valid);
    assert_abs_diff_eq!(results.coefficients.lift, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(results.coefficients.side, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(results.far_field_lift, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(results.far_field_induced_drag, 0.0, epsilon = TOL);
}

/// Ports `PositiveAlphaGivesPositiveLift`.
#[test]
fn positive_alpha_gives_positive_lift_and_drag() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let (_, results) = solve(&aircraft, 5.0);
    assert!(results.is_valid);
    assert!(results.coefficients.lift > 0.0);
    assert!(results.far_field_lift > 0.0);
    assert!(results.far_field_induced_drag > 0.0);
    assert!(results.span_efficiency > 0.0);
}

/// Ports `NegativeAlphaGivesNegativeLift`.
///
/// Induced drag is quadratic in lift, so it stays positive when the lift
/// changes sign — that asymmetry is the point of the case.
#[test]
fn negative_alpha_reverses_lift_but_not_induced_drag() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let (_, results) = solve(&aircraft, -5.0);
    assert!(results.is_valid);
    assert!(results.coefficients.lift < 0.0);
    assert!(results.far_field_lift < 0.0);
    assert!(results.far_field_induced_drag > 0.0);
}

/// Ports `LiftCoefficientMatchesCirculation`.
///
/// The far-field lift is recomputed here straight from the strip circulations
/// via the Kutta-Joukowski sum, independently of the Trefftz-plane routine, so
/// agreement is a real cross-check rather than a restatement.
#[test]
fn far_field_lift_matches_the_circulation_sum() {
    let span = 4.0;
    let chord = 1.0;
    let aircraft = rectangular_wing(span, chord, 4, 8);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();

    let reference_area = aircraft.reference_area();
    let mut lift_from_circulation = 0.0;
    for strip in solver.mesh().strips() {
        let circulation: f64 = (0..strip.num_vortices)
            .map(|k| solution.circulation[strip.first_vortex_index + k])
            .sum();
        lift_from_circulation += 2.0 * circulation * strip.width / reference_area;
    }

    let results = ForceCalculator::new().compute_for_aircraft(
        &aircraft,
        solver.mesh(),
        &solution,
        &condition,
    );
    assert_abs_diff_eq!(
        results.far_field_lift,
        lift_from_circulation,
        epsilon = 0.01 * lift_from_circulation.abs()
    );
}

/// Ports `SpanEfficiencyReasonable`.
#[test]
fn span_efficiency_of_a_rectangular_wing_is_near_elliptic() {
    let aircraft = rectangular_wing(10.0, 1.0, 4, 20);
    let (_, results) = solve(&aircraft, 5.0);
    assert!(results.is_valid);
    assert!(results.span_efficiency > 0.8);
    // Slightly above one is allowed: the discrete Trefftz integral overshoots
    // the ideal elliptic value a little at finite panel counts.
    assert!(results.span_efficiency < 1.15);
}

/// Ports `InducedDragQuadraticWithLift`.
///
/// Doubling alpha doubles lift and quadruples induced drag. This is the single
/// strongest statement that the far-field integral is the right functional
/// form and not merely positive.
#[test]
fn induced_drag_is_quadratic_in_lift() {
    let aircraft = rectangular_wing(10.0, 1.0, 4, 16);
    let (_, low) = solve(&aircraft, 2.0);
    let (_, high) = solve(&aircraft, 4.0);
    assert_abs_diff_eq!(
        high.far_field_lift / low.far_field_lift,
        2.0,
        epsilon = 0.05
    );
    assert_abs_diff_eq!(
        high.far_field_induced_drag / low.far_field_induced_drag,
        4.0,
        epsilon = 0.2
    );
}

/// Ports `SurfaceLoadsAddUpToTotal` and `StripLoadsHaveCorrectCount`.
#[test]
fn surface_loads_sum_to_the_total_and_carry_every_strip() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let (_, results) = solve(&aircraft, 5.0);
    assert!(results.is_valid);
    assert_eq!(results.surface_loads.len(), 1);
    assert_abs_diff_eq!(
        results.surface_loads[0].coefficients.lift,
        results.coefficients.lift,
        epsilon = 0.01
    );
    assert_abs_diff_eq!(
        results.surface_loads[0].coefficients.drag,
        results.coefficients.drag,
        epsilon = 0.01
    );
    assert_eq!(results.surface_loads[0].strip_loads.len(), 8);
}

/// Ports `PitchingMomentSign`.
///
/// With the moment reference moved aft of the quarter chord, a lifting wing
/// pitches nose up. This is the case that catches a sign error in the moment
/// arm, which the default reference at the origin would hide.
#[test]
fn an_aft_moment_reference_gives_positive_pitching_moment() {
    let mut aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    aircraft.set_moment_reference_point(Vec3::new(0.75, 0.0, 0.0));
    let (_, results) = solve(&aircraft, 5.0);
    assert!(results.is_valid);
    assert!(results.coefficients.pitch > 0.0);
}

/// Ports `AspectRatioCalculation`.
#[test]
fn results_report_the_geometric_aspect_ratio() {
    let aircraft = rectangular_wing(10.0, 1.0, 4, 8);
    let (_, results) = solve(&aircraft, 5.0);
    assert_abs_diff_eq!(results.aspect_ratio(), 10.0, epsilon = TOL);
}

/// Ports `IdealInducedDragComputation`.
#[test]
fn computed_induced_drag_brackets_the_elliptic_ideal() {
    let aircraft = rectangular_wing(10.0, 1.0, 4, 16);
    let (_, results) = solve(&aircraft, 5.0);
    let ideal = results.ideal_induced_drag();
    assert_abs_diff_eq!(
        ideal,
        results.coefficients.lift.powi(2) / (PI * results.aspect_ratio()),
        epsilon = TOL
    );
    // A rectangular planform is slightly worse than elliptic, never better by
    // more than numerical slack.
    assert!(results.far_field_induced_drag > ideal * 0.95);
    assert!(results.far_field_induced_drag < ideal * 1.25);
}

// ============================================================================
// stability_derivatives_test.cpp
// ============================================================================

fn derivatives(aircraft: &Aircraft, alpha_degrees: f64) -> zephyr_vlm::aero::StabilityDerivatives {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let condition = FlightCondition {
        alpha: alpha_degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    solver.solve(&condition).unwrap();
    StabilityDerivativeCalculator
        .compute(aircraft, &mut solver, &condition)
        .unwrap()
}

/// Ports `LiftCurveSlopePositive`.
#[test]
fn lift_curve_slope_is_positive() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let derivs = derivatives(&aircraft, 5.0);
    assert!(derivs.is_valid);
    assert!(derivs.by_alpha.lift > 0.0);
}

/// Ports `LiftCurveSlopeReasonable`.
///
/// Compared against the Helmbold/Prandtl finite-wing estimate
/// `2 pi / (1 + 2/AR)`, which a lifting-line-consistent VLM should land within
/// 30% of at AR 10.
#[test]
fn lift_curve_slope_matches_finite_wing_theory() {
    let (span, chord) = (10.0, 1.0);
    let aircraft = rectangular_wing(span, chord, 4, 16);
    let derivs = derivatives(&aircraft, 5.0);
    assert!(derivs.is_valid);
    let aspect_ratio = span / chord;
    let theory = 2.0 * PI / (1.0 + 2.0 / aspect_ratio);
    assert_abs_diff_eq!(derivs.by_alpha.lift, theory, epsilon = 0.3 * theory);
}

/// Ports `ZeroAlphaSymmetry`.
#[test]
fn a_symmetric_wing_has_no_side_force_response_to_alpha() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let derivs = derivatives(&aircraft, 0.0);
    assert!(derivs.is_valid);
    assert_abs_diff_eq!(derivs.by_alpha.side, 0.0, epsilon = TOL);
}

/// Ports `RollDampingNegative`, `PitchDampingFinite` and `YawDampingNonPositive`.
///
/// The damping signs are what make the aircraft dynamically stable; a positive
/// `Cl_p` would mean roll divergence.
#[test]
fn damping_derivatives_have_stable_signs() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let derivs = derivatives(&aircraft, 5.0);
    assert!(derivs.is_valid);
    assert!(derivs.by_roll_rate.roll < 0.0);
    assert!(derivs.by_pitch_rate.pitch.is_finite());
    assert!(derivs.by_yaw_rate.yaw <= TOL);
}

/// Ports `BodyAxisDerivativesValid`.
#[test]
fn body_axis_derivatives_are_populated() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    solver.solve(&condition).unwrap();
    let derivs = StabilityDerivativeCalculator
        .compute_body_axis(&aircraft, &mut solver, &condition)
        .unwrap();
    assert!(derivs.is_valid);
    // CZ_alpha: a positive alpha increment loads the wing, and body-axis Z
    // points down, so the response is positive in this convention.
    assert!(derivs.by_alpha_force.z > 0.0);
}

/// Ports `DerivativesAtDifferentAlpha`.
///
/// A linear method's derivatives should be nearly independent of the operating
/// point; drifting more than 20% between 2 and 8 degrees would mean the unit
/// responses are being contaminated by the base solution.
#[test]
fn derivatives_are_nearly_independent_of_alpha() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let low = derivatives(&aircraft, 2.0);
    let high = derivatives(&aircraft, 8.0);
    assert_abs_diff_eq!(
        low.by_alpha.lift,
        high.by_alpha.lift,
        epsilon = 0.2 * low.by_alpha.lift.abs()
    );
    assert_abs_diff_eq!(
        low.by_roll_rate.roll,
        high.by_roll_rate.roll,
        epsilon = 0.2 * low.by_roll_rate.roll.abs()
    );
}

/// Ports `StaticMarginComputed`.
#[test]
fn static_margin_is_finite() {
    let aircraft = rectangular_wing(4.0, 1.0, 4, 8);
    let derivs = derivatives(&aircraft, 5.0);
    assert!(derivs.is_valid);
    assert!(derivs.static_margin.is_finite());
}

/// Ports `ReferenceQuantitiesPreserved`.
#[test]
fn derivatives_carry_the_aircraft_reference_quantities() {
    let (span, chord) = (8.0, 2.0);
    let aircraft = rectangular_wing(span, chord, 4, 8);
    let derivs = derivatives(&aircraft, 5.0);
    assert_abs_diff_eq!(derivs.reference_area, span * chord, epsilon = TOL);
    assert_abs_diff_eq!(derivs.reference_chord, chord, epsilon = TOL);
    assert_abs_diff_eq!(derivs.reference_span, span, epsilon = TOL);
}

/// Ports `HighAspectRatioRollDamping`.
///
/// Roll damping grows with span, so the higher aspect ratio wing must be
/// *more* negative — a magnitude comparison that a sign-only test would miss.
#[test]
fn roll_damping_strengthens_with_aspect_ratio() {
    let low = derivatives(&rectangular_wing(4.0, 1.0, 4, 8), 5.0);
    let high = derivatives(&rectangular_wing(16.0, 1.0, 4, 32), 5.0);
    assert!(high.by_roll_rate.roll < low.by_roll_rate.roll);
}
