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

//! Newton convergence of the aircraft trim solver.
//!
//! The trim Jacobian once assembled each equation only against variables of
//! its own kind, dropping `dCL/d(elevator)` and `dCm/d(alpha)` from a
//! lift-and-elevator trim, and took its entries from the reported stability
//! derivatives, which hold the induced velocity frozen. The iteration still
//! reached the right answer but only linearly, about 0.3 per step, taking 12
//! iterations where Newton needs 3. These tests pin the exact Newton step
//! against finite differences of the residual the solver drives to zero,
//! the quadratic convergence that follows from it, and the converged state the
//! old iteration reached.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::{Coefficients, ForceCalculator, StabilityDerivativeCalculator};
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{Aircraft, ControlSurfaceDefinition, Section, Surface};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};
use zephyr_vlm::trim::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimOptions, TrimPresets, TrimResult,
    TrimSolver, TrimStatus, VariableConstraint,
};

/// Swept, dihedral, twisted wing with a tail carrying a full-span elevator:
/// 144 wing and 48 tail panels.
fn wing_and_tail() -> Aircraft {
    let elevator = ControlSurfaceDefinition {
        name: "elevator".into(),
        hinge_x_fraction: 0.7,
        hinge_vector: Vec3::y(),
        gain: 1.0,
        sign: 1.0,
    };
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(6);
    wing.set_spanwise_panel_count(12);
    wing.set_y_duplicate(0.0);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.5, 2.0));
    wing.add_section(Section::new(Vec3::new(0.3, 5.0, 0.2), 1.0, 0.0));
    let mut tail = Surface::new("Tail");
    tail.set_chordwise_panel_count(6);
    tail.set_spanwise_panel_count(4);
    tail.set_y_duplicate(0.0);
    let mut root = Section::new(Vec3::new(4.0, 0.0, 0.3), 0.6, 0.0);
    root.add_control_surface(elevator.clone());
    let mut tip = Section::new(Vec3::new(4.2, 1.5, 0.3), 0.4, 0.0);
    tip.add_control_surface(elevator);
    tail.add_section(root);
    tail.add_section(tip);
    let mut aircraft = Aircraft::new("wing-tail");
    aircraft.set_reference_area(12.5);
    aircraft.set_reference_chord(1.25);
    aircraft.set_reference_span(10.0);
    aircraft.set_moment_reference_point(Vec3::new(0.4, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft.add_surface(tail);
    aircraft
}

fn drives(kind: ConstraintType, target: f64) -> VariableConstraint {
    VariableConstraint { kind, target }
}

fn trim(
    aircraft: &Aircraft,
    constraints: &TrimConstraints,
    options: TrimOptions,
    initial: Option<&OperatingVariables>,
) -> TrimResult {
    let mut trim = TrimSolver::new();
    trim.set_options(options);
    let mut solver = VortexLatticeSolver::new();
    trim.solve(aircraft, &mut solver, constraints, initial)
}

fn elevator(state: &OperatingVariables) -> f64 {
    state
        .control_deflections
        .iter()
        .find(|(name, _)| name == "elevator")
        .map(|(_, value)| *value)
        .expect("elevator deflection")
}

/// alpha, beta, p, q, r, elevator.
fn variables(state: &OperatingVariables) -> [f64; 6] {
    [
        state.alpha,
        state.beta,
        state.roll_rate,
        state.pitch_rate,
        state.yaw_rate,
        elevator(state),
    ]
}

/// The reported coefficients the trim residual is built from.
fn coefficients(
    aircraft: &Aircraft,
    solver: &mut VortexLatticeSolver,
    x: [f64; 6],
) -> Coefficients {
    let condition = FlightCondition {
        alpha: x[0],
        beta: x[1],
        roll_rate: x[2],
        pitch_rate: x[3],
        yaw_rate: x[4],
        control_deflections: vec![("elevator".into(), x[5])],
        ..FlightCondition::default()
    };
    let solution = solver.solve(&condition).unwrap();
    ForceCalculator::default()
        .compute_for_aircraft(aircraft, solver.mesh(), &solution, &condition)
        .coefficients
}

/// The first Newton step solves `J delta = -r`. Checking it against a
/// central-difference Jacobian of the actual residual pins every column the
/// step moves along: an entry `J[i][j]` in error by `e` leaves row `i` off by
/// `e * delta[j]`.
///
/// All six variables are active and every equation couples to every variable
/// here (sideslip, rates, and an elevator at once), so the state rows meet the
/// control column and the control row meets the state columns.
#[test]
fn the_newton_step_solves_the_finite_difference_linearization() {
    let aircraft = wing_and_tail();
    let constraints = TrimConstraints {
        alpha: drives(ConstraintType::Lift, 0.5),
        beta: drives(ConstraintType::Side, 0.0),
        roll_rate: drives(ConstraintType::RollMoment, 0.0),
        pitch_rate: drives(ConstraintType::Fixed, 0.0),
        yaw_rate: drives(ConstraintType::YawMoment, 0.0),
        controls: vec![("elevator".into(), drives(ConstraintType::PitchMoment, 0.0))],
    };
    let initial = OperatingVariables {
        alpha: 0.06,
        beta: 0.02,
        roll_rate: 0.01,
        pitch_rate: 0.02,
        yaw_rate: 0.01,
        control_deflections: vec![("elevator".into(), 0.03)],
    };
    let options = TrimOptions {
        max_iterations: 1,
        tolerance: 0.0,
        // Only the step is under test; an inexact Jacobian may take a long
        // one and should fail on the rows below, not on a bound.
        max_alpha: f64::INFINITY,
        max_beta: f64::INFINITY,
        max_rate: f64::INFINITY,
        ..TrimOptions::default()
    };
    let stepped = trim(&aircraft, &constraints, options, Some(&initial));
    assert_eq!(
        stepped.status,
        TrimStatus::MaximumIterations,
        "{}",
        stepped.diagnostic
    );

    let x0 = variables(&initial);
    let x1 = variables(&stepped.state);
    let delta: Vec<f64> = x1
        .iter()
        .zip(x0)
        .map(|(after, before)| after - before)
        .collect();
    // The fixed pitch-rate row is satisfied exactly by the step.
    assert_abs_diff_eq!(x1[3], 0.0, epsilon = 1.0e-15);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    // CL, CY, Cl, Cn, Cm: the five coefficient equations.
    let outputs = |c: Coefficients| [c.lift, c.side, c.roll, c.yaw, c.pitch];
    let targets = [0.5, 0.0, 0.0, 0.0, 0.0];
    let base = outputs(coefficients(&aircraft, &mut solver, x0));
    let residual: Vec<f64> = base.iter().zip(targets).map(|(c, t)| c - t).collect();
    let step = 1.0e-6;
    let mut jacobian = [[0.0; 6]; 5];
    for (column, _) in x0.iter().enumerate() {
        let mut plus = x0;
        plus[column] += step;
        let mut minus = x0;
        minus[column] -= step;
        let up = outputs(coefficients(&aircraft, &mut solver, plus));
        let down = outputs(coefficients(&aircraft, &mut solver, minus));
        for row in 0..5 {
            jacobian[row][column] = (up[row] - down[row]) / (2.0 * step);
        }
    }
    let scale = residual.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    assert!(scale > 0.1, "the initial state should be far from trim");
    for row in 0..5 {
        let linearized: f64 = (0..6).map(|j| jacobian[row][j] * delta[j]).sum();
        // Central differences at this step are accurate to about 1e-9; the
        // stability-derivative Jacobian missed by 1e-3 of the residual and the
        // block-diagonal one by far more.
        assert!(
            (linearized + residual[row]).abs() <= 1.0e-7 * scale,
            "row {row}: J_fd delta = {linearized:e}, -r = {:e}",
            -residual[row]
        );
    }
}

/// Converged state of the pre-fix iteration on [`wing_and_tail`] at the
/// default tolerance, and with the tolerance tightened to 1e-14 (26
/// iterations).
const OLD_DEFAULT_ALPHA: f64 = 8.605_821_889_678_428e-2;
const OLD_DEFAULT_ELEVATOR: f64 = -7.772_497_685_799_65e-2;
const OLD_TIGHT_ALPHA: f64 = 8.605_824_400_786_438e-2;
const OLD_TIGHT_ELEVATOR: f64 = -7.772_530_190_630_907e-2;

#[test]
fn an_elevator_trim_converges_in_a_few_newton_iterations() {
    let aircraft = wing_and_tail();
    let constraints = TrimPresets::level_trimmed_flight(0.5, "elevator");
    let options = TrimOptions::default();
    let result = trim(&aircraft, &constraints, options, None);
    assert!(result.converged, "{}", result.diagnostic);
    // Newton reaches the default tolerance in 3; the old iteration took 12.
    assert!(result.iterations <= 6, "{} iterations", result.iterations);
    assert!(result.residual_norm < options.tolerance);
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.lift,
        0.5,
        epsilon = 1.0e-6
    );
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.pitch,
        0.0,
        epsilon = 1.0e-6
    );
    // Both iterations stop within the tolerance of the same root.
    assert_abs_diff_eq!(
        result.state.alpha,
        OLD_DEFAULT_ALPHA,
        epsilon = options.tolerance
    );
    assert_abs_diff_eq!(
        elevator(&result.state),
        OLD_DEFAULT_ELEVATOR,
        epsilon = options.tolerance
    );
}

#[test]
fn a_tight_elevator_trim_reaches_the_same_root_as_before() {
    let aircraft = wing_and_tail();
    let constraints = TrimPresets::level_trimmed_flight(0.5, "elevator");
    let options = TrimOptions {
        tolerance: 1.0e-13,
        ..TrimOptions::default()
    };
    let result = trim(&aircraft, &constraints, options, None);
    assert!(result.converged, "{}", result.diagnostic);
    assert!(result.iterations <= 6, "{} iterations", result.iterations);
    assert_abs_diff_eq!(result.state.alpha, OLD_TIGHT_ALPHA, epsilon = 1.0e-12);
    assert_abs_diff_eq!(
        elevator(&result.state),
        OLD_TIGHT_ELEVATOR,
        epsilon = 1.0e-12
    );
}

/// Residual after each of the first few iterations: the residual reported on
/// stopping at the iteration limit is the one evaluated at the last state.
fn residual_history(constraints: &TrimConstraints, count: usize) -> Vec<f64> {
    let aircraft = wing_and_tail();
    (1..=count)
        .map(|max_iterations| {
            let options = TrimOptions {
                max_iterations,
                tolerance: 0.0,
                ..TrimOptions::default()
            };
            trim(&aircraft, constraints, options, None).residual_norm
        })
        .collect()
}

/// Quadratic convergence, `r[k+1] <= C r[k]^2`. The measured constants are
/// about 0.015 and 0.07 for the two steps here; a linear iteration with even a
/// 0.1 per cent Jacobian error fails the second step by orders of magnitude.
#[test]
fn trim_residuals_fall_quadratically() {
    for constraints in [
        TrimPresets::level_trimmed_flight(0.5, "elevator"),
        TrimPresets::level_flight(0.5),
    ] {
        let history = residual_history(&constraints, 3);
        assert!(history[0] > 0.1);
        for pair in history.windows(2) {
            assert!(
                pair[1] <= pair[0] * pair[0],
                "residual history {history:?} is not quadratic"
            );
        }
    }
}

/// The trim reports derivatives in the same frozen-induced-velocity
/// convention as before, evaluated at the converged state: the Newton
/// Jacobian no longer comes from them, but the reported values must not move.
#[test]
fn the_trim_reports_the_stability_derivatives_of_its_converged_state() {
    let aircraft = wing_and_tail();
    let result = trim(
        &aircraft,
        &TrimPresets::level_trimmed_flight(0.5, "elevator"),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    let condition = FlightCondition {
        alpha: result.state.alpha,
        control_deflections: result.state.control_deflections.clone(),
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let expected = StabilityDerivativeCalculator
        .compute(&aircraft, &mut solver, &condition)
        .unwrap();
    let reported = &result.derivatives;
    for (actual, wanted) in [
        (reported.by_alpha.lift, expected.by_alpha.lift),
        (reported.by_alpha.pitch, expected.by_alpha.pitch),
        (reported.by_pitch_rate.pitch, expected.by_pitch_rate.pitch),
        (
            reported.control_derivatives[0].coefficients.pitch,
            expected.control_derivatives[0].coefficients.pitch,
        ),
    ] {
        assert_abs_diff_eq!(actual, wanted, epsilon = 1.0e-12);
    }
}

/// A constrained control the aircraft does not declare has no derivative. The
/// step is still taken, and the next aerodynamic solve names the control.
#[test]
fn an_undeclared_trim_control_is_reported_by_the_aerodynamic_solve() {
    let aircraft = wing_and_tail();
    let result = trim(
        &aircraft,
        &TrimPresets::level_trimmed_flight(0.5, "aileron"),
        TrimOptions::default(),
        None,
    );
    assert_eq!(result.status, TrimStatus::AerodynamicFailure);
    assert_eq!(result.iterations, 2);
    assert!(
        result.diagnostic.contains("unknown control 'aileron'"),
        "{}",
        result.diagnostic
    );
}

/// Derivatives are computed from a solution the solver already holds only
/// when it is still current; after the reference point moves they must match
/// a fresh solver's.
#[test]
fn stability_derivatives_resolve_after_the_reference_point_moves() {
    let aircraft = wing_and_tail();
    let condition = FlightCondition {
        alpha: 0.05,
        pitch_rate: 0.01,
        control_deflections: vec![("elevator".into(), 0.02)],
        ..FlightCondition::default()
    };
    let point = Vec3::new(1.0, 0.0, 0.1);

    let mut reused = VortexLatticeSolver::new();
    reused.set_aircraft(&aircraft).unwrap();
    reused.solve(&condition).unwrap();
    let first = StabilityDerivativeCalculator
        .compute(&aircraft, &mut reused, &condition)
        .unwrap();
    let again = StabilityDerivativeCalculator
        .compute(&aircraft, &mut reused, &condition)
        .unwrap();
    assert_eq!(first.by_alpha.lift, again.by_alpha.lift);
    assert_eq!(first.by_pitch_rate.pitch, again.by_pitch_rate.pitch);
    reused.set_reference_point(point);
    let moved = StabilityDerivativeCalculator
        .compute(&aircraft, &mut reused, &condition)
        .unwrap();
    let body = StabilityDerivativeCalculator
        .compute_body_axis(&aircraft, &mut reused, &condition)
        .unwrap();

    let mut fresh = VortexLatticeSolver::new();
    fresh.set_aircraft(&aircraft).unwrap();
    fresh.set_reference_point(point);
    let expected = StabilityDerivativeCalculator
        .compute(&aircraft, &mut fresh, &condition)
        .unwrap();
    let expected_body = StabilityDerivativeCalculator
        .compute_body_axis(&aircraft, &mut fresh, &condition)
        .unwrap();
    assert_abs_diff_eq!(
        moved.by_pitch_rate.lift,
        expected.by_pitch_rate.lift,
        epsilon = 1.0e-12
    );
    assert_abs_diff_eq!(
        moved.by_pitch_rate.pitch,
        expected.by_pitch_rate.pitch,
        epsilon = 1.0e-12
    );
    assert_abs_diff_eq!(
        body.by_pitch_rate_moment.y,
        expected_body.by_pitch_rate_moment.y,
        epsilon = 1.0e-12
    );
}
