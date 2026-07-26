mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::DEGREES_TO_RADIANS;
use pegasus_vlm::solver::VortexLatticeSolver;
use pegasus_vlm::trim::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimOptions, TrimPresets, TrimSolver,
    TrimStatus, VariableConstraint,
};

use common::{controlled_wing, rectangular_wing};

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
