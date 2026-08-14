//! Ported from the frozen v0.2.0 `flight_constraints_test.cpp`.
//!
//! The counting rules here are what size the trim Newton system: a `Free`
//! variable contributes neither an unknown nor an equation, so
//! `variable_count` and `constraint_count` have to agree with each other or
//! the Jacobian comes out non-square.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::DEGREES_TO_RADIANS;
use pegasus_vlm::trim::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimPresets, VariableConstraint,
};

const TOL: f64 = 1.0e-10;

fn constrained(kind: ConstraintType, target: f64) -> VariableConstraint {
    VariableConstraint { kind, target }
}

/// Ports `OperatingVariablesTest.DefaultNumVariables` and
/// `NumVariablesWithControls`.
#[test]
fn operating_variable_count_tracks_the_control_list() {
    let mut variables = OperatingVariables::default();
    // alpha, beta and the three body rates.
    assert_eq!(variables.variable_count(), 5);
    variables.control_deflections.push(("aileron".into(), 0.0));
    assert_eq!(variables.variable_count(), 6);
    variables.control_deflections.push(("elevator".into(), 0.1));
    assert_eq!(variables.variable_count(), 7);
}

/// Ports `OperatingVariablesTest.DefaultValues`.
#[test]
fn operating_variables_default_to_the_origin() {
    let variables = OperatingVariables::default();
    assert_abs_diff_eq!(variables.alpha, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(variables.beta, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(variables.roll_rate, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(variables.pitch_rate, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(variables.yaw_rate, 0.0, epsilon = TOL);
    assert!(variables.control_deflections.is_empty());
}

/// Ports `TrimConstraintsTest.DefaultAllFree`.
#[test]
fn default_constraints_are_all_free_and_count_zero() {
    let constraints = TrimConstraints::default();
    assert_eq!(constraints.alpha.kind, ConstraintType::Free);
    assert_eq!(constraints.beta.kind, ConstraintType::Free);
    assert_eq!(constraints.roll_rate.kind, ConstraintType::Free);
    assert_eq!(constraints.pitch_rate.kind, ConstraintType::Free);
    assert_eq!(constraints.yaw_rate.kind, ConstraintType::Free);
    assert_eq!(constraints.constraint_count(), 0);
}

/// Ports `TrimConstraintsTest.ConstraintCountWithControls`.
///
/// A `Free` control is listed but not solved for, so appending one must leave
/// the count alone — that is the case a naive `controls.len()` gets wrong.
#[test]
fn constraint_count_ignores_free_controls() {
    let mut constraints = TrimConstraints {
        alpha: constrained(ConstraintType::Lift, 0.5),
        ..TrimConstraints::default()
    };
    constraints.controls.push((
        "elevator".into(),
        constrained(ConstraintType::PitchMoment, 0.0),
    ));
    assert_eq!(constraints.constraint_count(), 2);

    constraints
        .controls
        .push(("rudder".into(), VariableConstraint::default()));
    assert_eq!(constraints.constraint_count(), 2);
}

/// Ports `TrimConstraintsTest.AllConstraintTypesCount`.
#[test]
fn every_state_constraint_type_counts_once() {
    let constraints = TrimConstraints {
        alpha: constrained(ConstraintType::Lift, 0.0),
        beta: constrained(ConstraintType::Side, 0.0),
        roll_rate: constrained(ConstraintType::RollMoment, 0.0),
        pitch_rate: constrained(ConstraintType::PitchMoment, 0.0),
        yaw_rate: constrained(ConstraintType::YawMoment, 0.0),
        controls: Vec::new(),
    };
    assert_eq!(constraints.constraint_count(), 5);
}

/// Ports `TrimPresetsTest.LevelFlightRatesFixed`.
#[test]
fn level_flight_preset_pins_sideslip_and_all_three_rates() {
    let constraints = TrimPresets::level_flight(0.5);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.5, epsilon = TOL);
    for pinned in [
        &constraints.beta,
        &constraints.roll_rate,
        &constraints.pitch_rate,
        &constraints.yaw_rate,
    ] {
        assert_eq!(pinned.kind, ConstraintType::Fixed);
        assert_abs_diff_eq!(pinned.target, 0.0, epsilon = TOL);
    }
    assert_eq!(constraints.constraint_count(), 5);
}

/// Ports `TrimPresetsTest.LevelTrimmedFlight`.
#[test]
fn level_trimmed_flight_preset_adds_the_named_control() {
    let constraints = TrimPresets::level_trimmed_flight(0.4, "elevator");
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.4, epsilon = TOL);
    assert_eq!(constraints.controls.len(), 1);
    assert_eq!(constraints.controls[0].0, "elevator");
    assert_eq!(constraints.controls[0].1.kind, ConstraintType::PitchMoment);
    assert_abs_diff_eq!(constraints.controls[0].1.target, 0.0, epsilon = TOL);
    assert_eq!(constraints.constraint_count(), 6);
}

/// Ports `TrimPresetsTest.SteadyTurn`.
///
/// The turn preset deliberately leaves pitch and yaw rate free — they are the
/// unknowns the turn solves for — so it carries three constraints, not five.
#[test]
fn steady_turn_preset_leaves_pitch_and_yaw_rate_free() {
    let constraints = TrimPresets::steady_turn(0.6, 30.0 * DEGREES_TO_RADIANS);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.6, epsilon = TOL);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.beta.target, 0.0, epsilon = TOL);
    assert_eq!(constraints.roll_rate.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.roll_rate.target, 0.0, epsilon = TOL);
    assert_eq!(constraints.pitch_rate.kind, ConstraintType::Free);
    assert_eq!(constraints.yaw_rate.kind, ConstraintType::Free);
    assert_eq!(constraints.constraint_count(), 3);
}

/// Ports `TrimPresetsTest.SteadySideslipAllFields`.
#[test]
fn steady_sideslip_preset_holds_the_requested_sideslip() {
    let sideslip = 5.0 * DEGREES_TO_RADIANS;
    let constraints = TrimPresets::steady_sideslip(0.3, sideslip);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.3, epsilon = TOL);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.beta.target, sideslip, epsilon = TOL);
    assert_eq!(constraints.roll_rate.kind, ConstraintType::Fixed);
    assert_eq!(constraints.pitch_rate.kind, ConstraintType::Fixed);
    assert_eq!(constraints.yaw_rate.kind, ConstraintType::Fixed);
    assert_eq!(constraints.constraint_count(), 5);
}

/// Ports `VariableConstraintTest.DefaultState`.
#[test]
fn variable_constraint_defaults_to_free_at_zero() {
    let constraint = VariableConstraint::default();
    assert_eq!(constraint.kind, ConstraintType::Free);
    assert_abs_diff_eq!(constraint.target, 0.0, epsilon = TOL);
}
