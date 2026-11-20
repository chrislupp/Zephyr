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

mod common;

use approx::assert_abs_diff_eq;
use zephyr_vlm::core::DEGREES_TO_RADIANS;
use zephyr_vlm::solver::VortexLatticeSolver;
use zephyr_vlm::trim::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimOptions, TrimPresets, TrimSolver,
    TrimStatus, VariableConstraint,
};

use common::{controlled_wing, lateral_directional_aircraft, rectangular_wing};

fn fixed(target: f64) -> VariableConstraint {
    VariableConstraint {
        kind: ConstraintType::Fixed,
        target,
    }
}

fn drives(kind: ConstraintType, target: f64) -> VariableConstraint {
    VariableConstraint { kind, target }
}

fn coefficient(results: &zephyr_vlm::aero::AerodynamicResults, kind: ConstraintType) -> f64 {
    match kind {
        ConstraintType::Lift => results.coefficients.lift,
        ConstraintType::Side => results.coefficients.side,
        ConstraintType::RollMoment => results.coefficients.roll,
        ConstraintType::PitchMoment => results.coefficients.pitch,
        ConstraintType::YawMoment => results.coefficients.yaw,
        ConstraintType::Free | ConstraintType::Fixed => {
            panic!("{kind:?} does not name a coefficient")
        }
    }
}

fn deflection(state: &OperatingVariables, name: &str) -> Option<f64> {
    state
        .control_deflections
        .iter()
        .find(|(candidate, _)| candidate == name)
        .map(|(_, value)| *value)
}

#[test]
fn operating_variables_count_includes_ordered_controls() {
    let state = OperatingVariables {
        control_deflections: vec![("a".into(), 0.0), ("a".into(), 0.1)],
        ..OperatingVariables::default()
    };
    assert_eq!(state.variable_count(), 7);
}

#[test]
fn trim_constraint_count_includes_each_active_control() {
    let constraints = TrimConstraints {
        alpha: VariableConstraint {
            kind: ConstraintType::Lift,
            target: 0.4,
        },
        beta: VariableConstraint {
            kind: ConstraintType::Fixed,
            target: 0.0,
        },
        controls: vec![
            (
                "elevator".into(),
                VariableConstraint {
                    kind: ConstraintType::PitchMoment,
                    target: 0.0,
                },
            ),
            ("free".into(), VariableConstraint::default()),
        ],
        ..TrimConstraints::default()
    };
    assert_eq!(constraints.constraint_count(), 3);
}

#[test]
fn level_flight_preset_fixes_beta_and_rates() {
    let constraints = TrimPresets::level_flight(0.5);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_eq!(constraints.alpha.target, 0.5);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_eq!(constraints.roll_rate.kind, ConstraintType::Fixed);
    assert_eq!(constraints.pitch_rate.kind, ConstraintType::Fixed);
    assert_eq!(constraints.yaw_rate.kind, ConstraintType::Fixed);
}

#[test]
fn steady_sideslip_preset_retains_requested_beta() {
    let constraints = TrimPresets::steady_sideslip(0.4, 3.0 * DEGREES_TO_RADIANS);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.beta.target, 3.0 * DEGREES_TO_RADIANS);
}

#[test]
fn trim_to_positive_lift_converges() {
    let aircraft = rectangular_wing(3, 8);
    let constraints = TrimPresets::level_flight(0.35);
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &constraints, None);
    assert!(result.converged, "{}", result.diagnostic);
    assert_eq!(result.status, TrimStatus::Converged);
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.lift,
        0.35,
        epsilon = 1.0e-6
    );
    assert!(result.state.alpha > 0.0);
    assert!(result.derivatives.is_valid);
}

#[test]
fn trim_to_zero_lift_converges_near_zero_alpha() {
    let aircraft = rectangular_wing(3, 6);
    let constraints = TrimPresets::level_flight(0.0);
    let initial = OperatingVariables {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..OperatingVariables::default()
    };
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &constraints, Some(&initial));
    assert!(result.converged);
    assert!(result.state.alpha.abs() < 1.0e-6);
}

#[test]
fn fixed_alpha_constraint_is_satisfied_exactly() {
    let aircraft = rectangular_wing(2, 4);
    let constraints = TrimConstraints {
        alpha: VariableConstraint {
            kind: ConstraintType::Fixed,
            target: 2.0 * DEGREES_TO_RADIANS,
        },
        ..TrimConstraints::default()
    };
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &constraints, None);
    assert!(result.converged);
    assert_abs_diff_eq!(
        result.state.alpha,
        2.0 * DEGREES_TO_RADIANS,
        epsilon = 1.0e-12
    );
}

#[test]
fn no_constraints_is_a_normal_nonconverged_result() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &TrimConstraints::default(), None);
    assert!(!result.converged);
    assert_eq!(result.status, TrimStatus::NoConstraints);
    assert!(result.diagnostic.contains("no constraints"));
}

#[test]
fn maximum_iteration_status_is_reported_without_panicking() {
    let aircraft = rectangular_wing(2, 4);
    let constraints = TrimPresets::level_flight(0.4);
    let mut trim = TrimSolver::new();
    trim.set_options(TrimOptions {
        max_iterations: 1,
        tolerance: 1.0e-16,
        relaxation: 0.1,
        ..TrimOptions::default()
    });
    let mut solver = VortexLatticeSolver::new();
    let result = trim.solve(&aircraft, &mut solver, &constraints, None);
    assert_eq!(result.status, TrimStatus::MaximumIterations);
    assert!(!result.converged);
    assert_eq!(result.iterations, 1);
}

#[test]
fn unreachable_lift_hits_configured_bounds() {
    let aircraft = rectangular_wing(2, 4);
    let constraints = TrimPresets::level_flight(10.0);
    let mut trim = TrimSolver::new();
    trim.set_options(TrimOptions {
        max_alpha: 0.1,
        ..TrimOptions::default()
    });
    let mut solver = VortexLatticeSolver::new();
    let result = trim.solve(&aircraft, &mut solver, &constraints, None);
    assert_eq!(result.status, TrimStatus::BoundsExceeded);
    assert!(!result.converged);
}

#[test]
fn controlled_pitch_trim_uses_control_derivative_column() {
    let aircraft = controlled_wing();
    let constraints = TrimConstraints {
        alpha: VariableConstraint {
            kind: ConstraintType::Fixed,
            target: 4.0 * DEGREES_TO_RADIANS,
        },
        controls: vec![(
            "elevator".into(),
            VariableConstraint {
                kind: ConstraintType::PitchMoment,
                target: 0.0,
            },
        )],
        ..TrimConstraints::default()
    };
    let initial = OperatingVariables {
        control_deflections: vec![("elevator".into(), 0.0)],
        ..OperatingVariables::default()
    };
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &constraints, Some(&initial));
    assert!(
        matches!(
            result.status,
            TrimStatus::Converged | TrimStatus::MaximumIterations
        ),
        "{}",
        result.diagnostic
    );
    assert!(
        result
            .state
            .control_deflections
            .iter()
            .any(|(name, value)| name == "elevator" && value.is_finite())
    );
}

/// Side force, roll moment and yaw moment as constraints on the *state*
/// variables.
///
/// Only `Lift` and `PitchMoment` had ever appeared in a test, so the three
/// lateral-directional arms of `build_system` were never entered. Driving all
/// three at once also proves the rows are wired to distinct coefficients: if
/// any two selected the same one the 3x3 block would be rank-deficient and the
/// solve would report `SingularJacobian` instead of converging.
#[test]
fn lateral_state_constraints_trim_their_own_coefficients() {
    let aircraft = lateral_directional_aircraft();
    let constraints = TrimConstraints {
        alpha: fixed(3.0 * DEGREES_TO_RADIANS),
        beta: drives(ConstraintType::Side, 0.0),
        roll_rate: drives(ConstraintType::RollMoment, 0.0),
        pitch_rate: fixed(0.0),
        yaw_rate: drives(ConstraintType::YawMoment, 0.0),
        controls: Vec::new(),
    };
    // Displaced away from the answer so the Newton step is actually taken;
    // starting at the solution would satisfy the assertions without ever
    // solving the lateral block.
    let initial = OperatingVariables {
        alpha: 3.0 * DEGREES_TO_RADIANS,
        beta: 4.0 * DEGREES_TO_RADIANS,
        roll_rate: 0.05,
        yaw_rate: 0.05,
        ..OperatingVariables::default()
    };
    // Tightened well below the 1e-6 default so "the coefficient equals its
    // target" is an assertion about the constraint wiring rather than about
    // the default stopping tolerance.
    let mut trim = TrimSolver::new();
    trim.set_options(TrimOptions {
        tolerance: 1.0e-12,
        max_iterations: 40,
        ..TrimOptions::default()
    });
    let mut solver = VortexLatticeSolver::new();
    let result = trim.solve(&aircraft, &mut solver, &constraints, Some(&initial));

    assert!(result.converged, "{}", result.diagnostic);
    assert_eq!(result.status, TrimStatus::Converged);
    assert!(result.iterations > 1, "converged without taking a step");
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.side,
        0.0,
        epsilon = 1.0e-9
    );
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.roll,
        0.0,
        epsilon = 1.0e-9
    );
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.yaw,
        0.0,
        epsilon = 1.0e-9
    );
    // CY = Cl = Cn = 0 has a unique root at beta = p = r = 0 for this
    // configuration, so recovering it confirms the residuals were the
    // coefficients themselves and not some other quantity that also vanishes.
    assert_abs_diff_eq!(result.state.beta, 0.0, epsilon = 1.0e-9);
    assert_abs_diff_eq!(result.state.roll_rate, 0.0, epsilon = 1.0e-9);
    assert_abs_diff_eq!(result.state.yaw_rate, 0.0, epsilon = 1.0e-9);
    assert_abs_diff_eq!(
        result.state.alpha,
        3.0 * DEGREES_TO_RADIANS,
        epsilon = 1.0e-12
    );
}

/// Every coefficient constraint applied to a *control* rather than a state.
///
/// This is the `kind =>` arm of the control loop, which reads both the target
/// coefficient and the matching control derivative. Only `PitchMoment` was
/// covered, so a mismatch between the two lookups -- residual on one
/// coefficient, Jacobian column from another -- was undetectable for the
/// remaining four.
#[test]
fn control_constraints_trim_each_coefficient_to_its_target() {
    for (control, kind, target) in [
        ("aileron", ConstraintType::Lift, 0.30),
        ("rudder", ConstraintType::Side, 0.02),
        ("aileron", ConstraintType::RollMoment, 0.05),
        ("rudder", ConstraintType::YawMoment, 0.03),
    ] {
        let aircraft = lateral_directional_aircraft();
        let other = if control == "aileron" {
            "rudder"
        } else {
            "aileron"
        };
        let constraints = TrimConstraints {
            alpha: fixed(3.0 * DEGREES_TO_RADIANS),
            beta: fixed(0.0),
            roll_rate: fixed(0.0),
            pitch_rate: fixed(0.0),
            yaw_rate: fixed(0.0),
            // The unconstrained control stays in the list so the inactive-control
            // skip is exercised in both the assembly and the state update.
            controls: vec![
                (control.into(), drives(kind, target)),
                (other.into(), VariableConstraint::default()),
            ],
        };
        // No deflections seeded: the update has to append the control to the
        // state rather than accumulate into an existing entry.
        let initial = OperatingVariables {
            alpha: 3.0 * DEGREES_TO_RADIANS,
            ..OperatingVariables::default()
        };
        let mut trim = TrimSolver::new();
        trim.set_options(TrimOptions {
            tolerance: 1.0e-12,
            max_iterations: 40,
            ..TrimOptions::default()
        });
        let mut solver = VortexLatticeSolver::new();
        let result = trim.solve(&aircraft, &mut solver, &constraints, Some(&initial));

        assert!(
            result.converged,
            "{control} / {kind:?}: {}",
            result.diagnostic
        );
        assert_abs_diff_eq!(
            coefficient(&result.aerodynamic_results, kind),
            target,
            epsilon = 1.0e-9
        );
        let moved = deflection(&result.state, control)
            .unwrap_or_else(|| panic!("{control} was never added to the trimmed state"));
        assert!(
            moved.abs() > 1.0e-6,
            "{control} / {kind:?}: target reached without deflecting the control"
        );
        assert_eq!(
            deflection(&result.state, other),
            None,
            "the free control must not be driven"
        );
    }
}

/// A control pinned by `Fixed`, the control-loop counterpart of
/// [`fixed_alpha_constraint_is_satisfied_exactly`].
#[test]
fn fixed_control_constraint_is_satisfied_exactly() {
    let aircraft = lateral_directional_aircraft();
    let constraints = TrimConstraints {
        alpha: fixed(2.0 * DEGREES_TO_RADIANS),
        beta: fixed(0.0),
        roll_rate: fixed(0.0),
        pitch_rate: fixed(0.0),
        yaw_rate: fixed(0.0),
        controls: vec![("rudder".into(), fixed(5.0 * DEGREES_TO_RADIANS))],
    };
    let mut solver = VortexLatticeSolver::new();
    let result = TrimSolver::new().solve(&aircraft, &mut solver, &constraints, None);

    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(
        deflection(&result.state, "rudder").expect("rudder missing from the trimmed state"),
        5.0 * DEGREES_TO_RADIANS,
        epsilon = 1.0e-12
    );
    // A deflected rudder on this fixture produces side force, so the fixed
    // deflection reached the solve rather than only the state vector.
    assert!(result.aerodynamic_results.coefficients.side.abs() > 1.0e-4);
}
