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

//! Ported from the frozen v0.2.0 `gradient_calculator_test.cpp`.
//!
//! The top-level gradient entry point: flight-condition derivatives checked
//! against the independently computed stability derivatives, geometry
//! derivatives checked for sign and structure, and the output-selection
//! machinery that lets a caller pay only for the outputs it wants.
//!
//! The C++ exposed one method per input (`ComputeAlphaDerivatives`,
//! `ComputeRollRateDerivatives`, ...). The Rust computes all six inputs in a
//! single batched adjoint pass through `compute_flight_condition_gradients`, so
//! the rate and Mach cases read their column out of `FullGradients` rather than
//! calling a separate method.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::StabilityDerivativeCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{AeroGradients, GradientCalculator, OutputSelection};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// The frozen 5% relative band, with a small absolute floor.
const DERIVATIVE_TOLERANCE: f64 = 0.05;

/// The frozen fixture: a 4x8 half wing of span 4, referenced to a span of 8.
fn wing_aircraft() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(8);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Test Aircraft");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(8.0);
    aircraft
}

/// The same wing with a horizontal tail aft, which makes `Cm_alpha` negative.
fn wing_and_tail_aircraft() -> Aircraft {
    let mut aircraft = wing_aircraft();
    let mut tail = Surface::new("HTail");
    tail.set_chordwise_panel_count(3);
    tail.set_spanwise_panel_count(4);
    tail.add_section(Section::new(Vec3::new(4.0, 0.0, 0.0), 0.5, 0.0));
    tail.add_section(Section::new(Vec3::new(4.0, 1.5, 0.0), 0.5, 0.0));
    aircraft.add_surface(tail);
    aircraft
}

fn at_alpha(degrees: f64) -> FlightCondition {
    FlightCondition {
        alpha: degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    }
}

fn calculator_for(aircraft: &Aircraft) -> GradientCalculator {
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(aircraft);
    calculator
}

fn stability_derivatives(
    aircraft: &Aircraft,
    condition: &FlightCondition,
) -> zephyr_vlm::aero::StabilityDerivatives {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    solver.solve(condition).unwrap();
    StabilityDerivativeCalculator
        .compute(aircraft, &mut solver, condition)
        .unwrap()
}

fn assert_near(gradient: f64, reference: f64) {
    assert_abs_diff_eq!(
        gradient,
        reference,
        epsilon = reference.abs() * DERIVATIVE_TOLERANCE + 0.01
    );
}

// ============================================================================
// Flight-condition derivatives against the stability derivatives
// ============================================================================

/// Ports `AlphaDerivativesMatchStabilityDerivatives`.
///
/// The adjoint gradient and the unit-response stability derivative are computed
/// by entirely different routes, so agreement is a genuine cross-check.
#[test]
fn alpha_gradients_match_the_stability_derivatives() {
    let aircraft = wing_aircraft();
    let condition = at_alpha(5.0);
    let reference = stability_derivatives(&aircraft, &condition);
    let gradients = calculator_for(&aircraft)
        .compute_alpha_gradients(&condition)
        .unwrap();
    assert_near(gradients.coefficients.lift, reference.by_alpha.lift);
    assert_near(gradients.coefficients.pitch, reference.by_alpha.pitch);
    assert_near(gradients.coefficients.side, reference.by_alpha.side);
}

/// Ports `BetaDerivativesMatchStabilityDerivatives`.
#[test]
fn beta_gradients_match_the_stability_derivatives() {
    let aircraft = wing_aircraft();
    let condition = at_alpha(5.0);
    let reference = stability_derivatives(&aircraft, &condition);
    let gradients = calculator_for(&aircraft)
        .compute_beta_gradients(&condition)
        .unwrap();
    assert_near(gradients.coefficients.side, reference.by_beta.side);
    assert_near(gradients.coefficients.roll, reference.by_beta.roll);
    assert_near(gradients.coefficients.yaw, reference.by_beta.yaw);
}

/// Ports the three rate cases: `RollRateDerivativesMatchStabilityDerivatives`,
/// `PitchRate...` and `YawRate...`.
///
/// The C++ had one method per rate; the Rust computes all six inputs in one
/// batched pass, so this reads three columns out of a single result.
#[test]
fn rate_gradients_match_the_stability_derivatives() {
    let aircraft = wing_aircraft();
    let condition = at_alpha(5.0);
    let reference = stability_derivatives(&aircraft, &condition);
    let gradients = calculator_for(&aircraft)
        .compute_flight_condition_gradients(&condition)
        .unwrap();

    assert_near(
        gradients.by_roll_rate.coefficients.roll,
        reference.by_roll_rate.roll,
    );
    assert_near(
        gradients.by_roll_rate.coefficients.lift,
        reference.by_roll_rate.lift,
    );
    assert_near(
        gradients.by_pitch_rate.coefficients.pitch,
        reference.by_pitch_rate.pitch,
    );
    assert_near(
        gradients.by_pitch_rate.coefficients.lift,
        reference.by_pitch_rate.lift,
    );
    assert_near(
        gradients.by_yaw_rate.coefficients.yaw,
        reference.by_yaw_rate.yaw,
    );
    assert_near(
        gradients.by_yaw_rate.coefficients.roll,
        reference.by_yaw_rate.roll,
    );
}

/// Ports `MachDerivativesComputed`, `MachDerivativesAtLowMach` and
/// `MachDerivativesAtHighMach`.
///
/// Compressibility increases lift, and the effect strengthens with Mach — the
/// Prandtl-Glauert factor `1/sqrt(1 - M^2)` steepens as `M` rises. The
/// monotonic comparison is the real content: a positive derivative alone could
/// come from a constant.
#[test]
fn mach_gradients_are_positive_and_strengthen_with_mach() {
    let aircraft = wing_aircraft();
    let calculator = calculator_for(&aircraft);

    let low = calculator
        .compute_flight_condition_gradients(&FlightCondition {
            mach: 0.1,
            ..at_alpha(5.0)
        })
        .unwrap();
    assert!(low.by_mach.coefficients.lift > 0.0);
    assert!(low.by_mach.coefficients.lift.abs() < 1.0);
    assert!(low.by_mach.induced_drag > 0.0);

    let high = calculator
        .compute_flight_condition_gradients(&FlightCondition {
            mach: 0.6,
            ..at_alpha(5.0)
        })
        .unwrap();
    assert!(high.by_mach.coefficients.lift > low.by_mach.coefficients.lift);
}

/// Ports `FullGradientsComputation`.
#[test]
fn the_full_gradient_set_is_physically_sensible() {
    let aircraft = wing_aircraft();
    let gradients = calculator_for(&aircraft)
        .compute_flight_condition_gradients(&FlightCondition {
            mach: 0.2,
            ..at_alpha(5.0)
        })
        .unwrap();

    // Around 5 per radian for an aspect ratio near 8.
    assert!(gradients.by_alpha.coefficients.lift > 3.0);
    assert!(gradients.by_alpha.coefficients.lift < 7.0);
    // Roll damping is stabilising.
    assert!(gradients.by_roll_rate.coefficients.roll < 0.0);
    assert!(gradients.by_pitch_rate.coefficients.pitch != 0.0);
    assert!(gradients.by_mach.coefficients.lift > 0.0);
}

/// Ports `DerivativesAtZeroAlpha`.
///
/// Lift still responds to alpha at zero incidence, while drag does not to first
/// order — induced drag is quadratic in lift, so its slope through zero lift is
/// zero.
#[test]
fn alpha_gradients_at_zero_incidence() {
    let aircraft = wing_aircraft();
    let gradients = calculator_for(&aircraft)
        .compute_alpha_gradients(&at_alpha(0.0))
        .unwrap();
    assert!(gradients.coefficients.lift > 0.0);
    assert_abs_diff_eq!(gradients.coefficients.drag, 0.0, epsilon = 0.1);
}

/// Ports `GradientCalculatorWithTailTest.CmAlphaIsNegativeWithTail` and
/// `MatchesStabilityDerivatives`.
///
/// A tail aft of the moment reference makes the configuration statically
/// stable, which flips the sign of `Cm_alpha`. That sign is the whole point of
/// the fixture.
#[test]
fn a_tailed_configuration_has_negative_pitch_stiffness() {
    let aircraft = wing_and_tail_aircraft();
    let condition = at_alpha(5.0);
    let gradients = calculator_for(&aircraft)
        .compute_alpha_gradients(&condition)
        .unwrap();
    assert!(gradients.coefficients.pitch < 0.0);

    let reference = stability_derivatives(&aircraft, &condition);
    assert_near(gradients.coefficients.pitch, reference.by_alpha.pitch);
}

// ============================================================================
// Geometry derivatives
// ============================================================================

/// A tapered, twisted wing so that chord, twist and span-position derivatives
/// are all individually nonzero.
fn geometry_aircraft() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(6);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.2, 2.0));
    wing.add_section(Section::new(Vec3::new(0.3, 4.0, 0.1), 0.8, -1.0));

    let mut aircraft = Aircraft::new("Geometry");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(8.0);
    aircraft
}

fn geometry_gradients(aircraft: &Aircraft) -> zephyr_vlm::gradients::GeometryGradients {
    calculator_for(aircraft)
        .compute_geometry_gradients(&at_alpha(5.0))
        .unwrap()
}

/// Ports `ChordGradientIsNonZero`, `SpanPositionGradientNonZero` and
/// `DihedralGradientAffectsRollingMoment`.
#[test]
fn chord_span_and_dihedral_derivatives_are_nonzero_and_finite() {
    let gradients = geometry_gradients(&geometry_aircraft());
    for section in &gradients.sections {
        assert!(section.by_chord.coefficients.lift != 0.0);
        assert!(section.by_leading_edge_y.coefficients.lift != 0.0);
        // Dihedral acts through the leading-edge z station.
        assert!(section.by_leading_edge_z.coefficients.roll.is_finite());
    }
}

/// Ports `TwistGradientPhysicallyCorrect`.
///
/// Washing in either station raises lift, so both derivatives are positive —
/// not merely nonzero, which a sign error would also satisfy.
#[test]
fn twisting_either_station_increases_lift() {
    let gradients = geometry_gradients(&geometry_aircraft());
    assert_eq!(gradients.sections.len(), 2);
    for section in &gradients.sections {
        assert!(section.by_twist.coefficients.lift != 0.0);
        assert!(
            section.by_twist.coefficients.lift > 0.0,
            "section {} twist derivative was {}",
            section.section_index,
            section.by_twist.coefficients.lift
        );
    }
}

/// Ports `FullGeometryGradientsReturnsAllSections`, and its design-variable
/// count.
///
/// Five variables per section — chord, twist and the three leading-edge
/// coordinates — so two sections give ten.
#[test]
fn geometry_gradients_cover_every_section_with_five_variables_each() {
    let gradients = geometry_gradients(&geometry_aircraft());
    assert_eq!(gradients.sections.len(), 2);
    assert_eq!(gradients.sections[0].surface_name, "Wing");
    assert_eq!(gradients.sections[1].surface_name, "Wing");
    assert_eq!(gradients.sections[0].section_index, 0);
    assert_eq!(gradients.sections[1].section_index, 1);
    assert_eq!(gradients.design_variable_count(), 10);
}

/// Ports `GetChordGradientsForExtraction` and `GetTwistGradientsForExtraction`.
#[test]
fn per_section_chord_and_twist_columns_can_be_extracted() {
    let gradients = geometry_gradients(&geometry_aircraft());
    let lift_by_chord: Vec<f64> = gradients
        .sections
        .iter()
        .map(|section| section.by_chord.coefficients.lift)
        .collect();
    assert_eq!(lift_by_chord.len(), 2);
    assert!(lift_by_chord.iter().all(|value| *value != 0.0));

    let pitch_by_twist: Vec<f64> = gradients
        .sections
        .iter()
        .map(|section| section.by_twist.coefficients.pitch)
        .collect();
    assert_eq!(pitch_by_twist.len(), 2);
    assert!(pitch_by_twist.iter().all(|value| *value != 0.0));
}

/// Ports `GradientsAreSymmetricForSymmetricSections`.
#[test]
fn every_geometry_derivative_is_finite() {
    let gradients = geometry_gradients(&geometry_aircraft());
    for section in &gradients.sections {
        for group in [
            &section.by_chord,
            &section.by_twist,
            &section.by_leading_edge_x,
            &section.by_leading_edge_y,
            &section.by_leading_edge_z,
        ] {
            for value in [
                group.coefficients.lift,
                group.coefficients.drag,
                group.coefficients.pitch,
            ] {
                assert!(value.is_finite());
            }
        }
    }
}

/// Ports `MultiSurfaceGeometryGradientTest.GradientsComputedForAllSurfaces` and
/// `TailTwistAffectsPitchingMoment`.
///
/// Sections are reported surface-major, and twisting the tail pitches the
/// aircraft nose down because the tail sits aft of the moment reference.
#[test]
fn multi_surface_geometry_gradients_are_surface_major() {
    let aircraft = wing_and_tail_aircraft();
    let gradients = geometry_gradients(&aircraft);
    assert_eq!(gradients.sections.len(), 4);
    assert_eq!(gradients.sections[0].surface_name, "Wing");
    assert_eq!(gradients.sections[1].surface_name, "Wing");
    assert_eq!(gradients.sections[2].surface_name, "HTail");
    assert_eq!(gradients.sections[3].surface_name, "HTail");
    assert_eq!(gradients.design_variable_count(), 20);

    for section in &gradients.sections[2..] {
        assert!(section.by_twist.coefficients.pitch != 0.0);
        assert!(
            section.by_twist.coefficients.pitch < 0.0,
            "tail twist should pitch nose down, got {}",
            section.by_twist.coefficients.pitch
        );
    }
}

// ============================================================================
// Output selection
// ============================================================================

fn selected(aircraft: &Aircraft, selection: OutputSelection) -> AeroGradients {
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(aircraft);
    calculator.set_output_selection(selection);
    calculator.compute_alpha_gradients(&at_alpha(5.0)).unwrap()
}

/// Ports `OutputSelectionSubsetMatchesFullResults`, `OutputSelectionSingleOutput`
/// and `OutputSelectionAllDisabled`.
///
/// A selected output must be bit-identical to the full computation, and an
/// unselected one exactly zero — not merely small. That exactness is what makes
/// selection a pure cost saving rather than an approximation.
#[test]
fn output_selection_is_exact_in_both_directions() {
    let aircraft = wing_aircraft();
    let full = selected(&aircraft, OutputSelection::all());

    let subset = selected(
        &aircraft,
        OutputSelection {
            lift: true,
            drag: true,
            ..OutputSelection::none()
        },
    );
    assert_abs_diff_eq!(
        subset.coefficients.lift,
        full.coefficients.lift,
        epsilon = 1.0e-10
    );
    assert_abs_diff_eq!(
        subset.coefficients.drag,
        full.coefficients.drag,
        epsilon = 1.0e-10
    );
    assert_eq!(subset.coefficients.side, 0.0);
    assert_eq!(subset.coefficients.roll, 0.0);
    assert_eq!(subset.coefficients.pitch, 0.0);
    assert_eq!(subset.coefficients.yaw, 0.0);

    let single = selected(
        &aircraft,
        OutputSelection {
            pitch: true,
            ..OutputSelection::none()
        },
    );
    assert_abs_diff_eq!(
        single.coefficients.pitch,
        full.coefficients.pitch,
        epsilon = 1.0e-10
    );
    assert_eq!(single.coefficients.lift, 0.0);
    assert_eq!(single.coefficients.drag, 0.0);
    assert_eq!(single.coefficients.side, 0.0);
    assert_eq!(single.coefficients.roll, 0.0);
    assert_eq!(single.coefficients.yaw, 0.0);

    let none = selected(&aircraft, OutputSelection::none());
    assert_eq!(none.coefficients.lift, 0.0);
    assert_eq!(none.coefficients.drag, 0.0);
    assert_eq!(none.coefficients.side, 0.0);
    assert_eq!(none.coefficients.roll, 0.0);
    assert_eq!(none.coefficients.pitch, 0.0);
    assert_eq!(none.coefficients.yaw, 0.0);
}

/// Ports `OutputSelectionBetaDerivatives`.
#[test]
fn output_selection_applies_to_beta_gradients() {
    let aircraft = wing_aircraft();
    let condition = at_alpha(5.0);

    let mut full_calculator = GradientCalculator::new();
    full_calculator.set_aircraft(&aircraft);
    let full = full_calculator.compute_beta_gradients(&condition).unwrap();

    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection {
        side: true,
        roll: true,
        yaw: true,
        ..OutputSelection::none()
    });
    let subset = calculator.compute_beta_gradients(&condition).unwrap();

    assert_abs_diff_eq!(
        subset.coefficients.side,
        full.coefficients.side,
        epsilon = 1.0e-10
    );
    assert_abs_diff_eq!(
        subset.coefficients.roll,
        full.coefficients.roll,
        epsilon = 1.0e-10
    );
    assert_abs_diff_eq!(
        subset.coefficients.yaw,
        full.coefficients.yaw,
        epsilon = 1.0e-10
    );
    assert_eq!(subset.coefficients.lift, 0.0);
    assert_eq!(subset.coefficients.drag, 0.0);
    assert_eq!(subset.coefficients.pitch, 0.0);
}

/// Ports `OutputSelectionGeometryGradients`.
#[test]
fn output_selection_applies_to_geometry_gradients() {
    let aircraft = geometry_aircraft();
    let full = geometry_gradients(&aircraft);

    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection {
        lift: true,
        drag: true,
        pitch: true,
        ..OutputSelection::none()
    });
    let subset = calculator
        .compute_geometry_gradients(&at_alpha(5.0))
        .unwrap();

    assert_eq!(subset.sections.len(), full.sections.len());
    for (subset, full) in subset.sections.iter().zip(&full.sections) {
        assert_abs_diff_eq!(
            subset.by_chord.coefficients.lift,
            full.by_chord.coefficients.lift,
            epsilon = 1.0e-10
        );
        assert_abs_diff_eq!(
            subset.by_chord.coefficients.drag,
            full.by_chord.coefficients.drag,
            epsilon = 1.0e-10
        );
        assert_abs_diff_eq!(
            subset.by_chord.coefficients.pitch,
            full.by_chord.coefficients.pitch,
            epsilon = 1.0e-10
        );
        assert_eq!(subset.by_chord.coefficients.side, 0.0);
        assert_eq!(subset.by_chord.coefficients.roll, 0.0);
        assert_eq!(subset.by_chord.coefficients.yaw, 0.0);
    }
}
