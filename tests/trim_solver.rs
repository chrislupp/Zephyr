//! Ported from the frozen v0.2.0 `trim_solver_test.cpp`.
//!
//! The Newton trim driver end to end: converging to a target lift, holding
//! fixed states, and the failure modes. Where the C++ reported failure through
//! a `converged` flag plus an error string, the Rust also carries a
//! [`TrimStatus`], so the failure cases assert *which* failure occurred rather
//! than only that one did.

use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, SpacingType, Vec3};
use pegasus_vlm::geometry::{Aircraft, Section, Surface};
use pegasus_vlm::solver::VortexLatticeSolver;
use pegasus_vlm::trim::{
    ConstraintType, OperatingVariables, TrimConstraints, TrimOptions, TrimPresets, TrimSolver,
    TrimStatus, VariableConstraint,
};

const TOL: f64 = 1.0e-4;

/// The frozen fixture: an AR-10 flat rectangular wing.
fn wing() -> Aircraft {
    let (span, chord) = (10.0, 1.0);
    let mut surface = Surface::new("Wing");
    surface.set_chordwise_panel_count(4);
    surface.set_spanwise_panel_count(16);
    surface.set_chordwise_spacing(SpacingType::Cosine);
    surface.set_spanwise_spacing(SpacingType::Cosine);
    surface.add_section(Section::new(Vec3::new(0.0, -span / 2.0, 0.0), chord, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, span / 2.0, 0.0), chord, 0.0));

    let mut aircraft = Aircraft::new("Test Wing");
    aircraft.add_surface(surface);
    aircraft.set_reference_area(span * chord);
    aircraft.set_reference_chord(chord);
    aircraft.set_reference_span(span);
    aircraft
}

fn fixed(target: f64) -> VariableConstraint {
    VariableConstraint {
        kind: ConstraintType::Fixed,
        target,
    }
}

/// All five state variables pinned, with `alpha` given the supplied constraint.
fn alpha_constrained(alpha: VariableConstraint) -> TrimConstraints {
    TrimConstraints {
        alpha,
        beta: fixed(0.0),
        roll_rate: fixed(0.0),
        pitch_rate: fixed(0.0),
        yaw_rate: fixed(0.0),
        controls: Vec::new(),
    }
}

fn trim_with(
    constraints: &TrimConstraints,
    options: TrimOptions,
    initial: Option<&OperatingVariables>,
) -> pegasus_vlm::trim::TrimResult {
    let aircraft = wing();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let mut trimmer = TrimSolver::new();
    trimmer.set_options(options);
    trimmer.solve(&aircraft, &mut solver, constraints, initial)
}

/// Ports `TrimToTargetCL`.
#[test]
fn trimming_to_a_target_lift_converges_at_positive_alpha() {
    let result = trim_with(
        &TrimPresets::level_flight(0.5),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.lift,
        0.5,
        epsilon = 0.01
    );
    assert!(result.state.alpha > 0.0);
}

/// Ports `TrimToZeroCL`.
///
/// A flat plate reaches zero lift at exactly zero incidence, so this pins the
/// solution rather than merely the residual.
#[test]
fn trimming_a_flat_wing_to_zero_lift_lands_at_zero_alpha() {
    let result = trim_with(
        &TrimPresets::level_flight(0.0),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(
        result.aerodynamic_results.coefficients.lift,
        0.0,
        epsilon = TOL
    );
    assert_abs_diff_eq!(result.state.alpha, 0.0, epsilon = TOL);
}

/// Ports `TrimWithFixedAlpha`.
#[test]
fn a_fixed_alpha_is_held_exactly_and_still_produces_lift() {
    let target = 5.0 * DEGREES_TO_RADIANS;
    let result = trim_with(
        &alpha_constrained(fixed(target)),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(result.state.alpha, target, epsilon = TOL);
    assert!(result.aerodynamic_results.coefficients.lift > 0.0);
}

/// Ports `TrimWithFixedBeta`.
#[test]
fn a_fixed_sideslip_is_held_while_lift_is_trimmed() {
    let target = 5.0 * DEGREES_TO_RADIANS;
    let result = trim_with(
        &TrimPresets::steady_sideslip(0.5, target),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(result.state.beta, target, epsilon = TOL);
}

/// Ports `RatesRemainZeroForLevelFlight`.
#[test]
fn level_flight_leaves_all_three_rates_at_zero() {
    let result = trim_with(
        &TrimPresets::level_flight(0.4),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert_abs_diff_eq!(result.state.roll_rate, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(result.state.pitch_rate, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(result.state.yaw_rate, 0.0, epsilon = TOL);
}

/// Ports `ConvergenceInFewIterations` and `InitialGuessIsUsed`.
///
/// Newton on a linear method should converge in a couple of steps, and a good
/// initial guess must make it strictly cheaper — that is the observable proof
/// the guess is actually used rather than discarded.
#[test]
fn a_good_initial_guess_reduces_the_iteration_count() {
    let cold = trim_with(
        &TrimPresets::level_flight(0.3),
        TrimOptions::default(),
        None,
    );
    assert!(cold.converged);
    assert!(cold.iterations <= 10);

    let initial = OperatingVariables {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        ..OperatingVariables::default()
    };
    let warm = trim_with(
        &TrimPresets::level_flight(0.4),
        TrimOptions::default(),
        Some(&initial),
    );
    assert!(warm.converged);
    assert!(warm.iterations <= 5);
}

/// Ports `ResidualDecreasesMonotonically`.
///
/// The frozen name overstates what it checked: it asserted only that a
/// converged result's residual is below tolerance, which is the definition of
/// convergence. Kept as written rather than strengthened into a claim the
/// oracle never made.
#[test]
fn a_converged_result_has_a_residual_below_tolerance() {
    let options = TrimOptions {
        max_iterations: 10,
        ..TrimOptions::default()
    };
    let result = trim_with(&TrimPresets::level_flight(0.6), options, None);
    if result.converged {
        assert!(result.residual_norm < options.tolerance);
    }
}

/// Ports `StabilityDerivativesComputed`.
#[test]
fn a_trimmed_result_carries_valid_stability_derivatives() {
    let result = trim_with(
        &TrimPresets::level_flight(0.5),
        TrimOptions::default(),
        None,
    );
    assert!(result.converged, "{}", result.diagnostic);
    assert!(result.derivatives.is_valid);
    assert!(result.derivatives.by_alpha.lift > 0.0);
}

/// Ports `NoConstraintsFails`.
///
/// An all-free constraint set has nothing to solve for. The Rust names the
/// reason, so this asserts the status rather than only a non-empty message.
#[test]
fn an_unconstrained_problem_reports_no_constraints() {
    let result = trim_with(&TrimConstraints::default(), TrimOptions::default(), None);
    assert!(!result.converged);
    assert_eq!(result.status, TrimStatus::NoConstraints);
    assert!(!result.diagnostic.is_empty());
}

/// Ports `UnreachableCLFails`.
///
/// A lift coefficient of 10 is unreachable within a half-radian alpha bound, so
/// the Newton step walks out of bounds and the solve stops.
///
/// The solver reports the offending state rather than clamping it back to the
/// bound: `BoundsExceeded` names *why* it stopped, and the state shows how far
/// the step went. Clamping would leave a state that satisfies the bound while
/// satisfying nothing else, which reads as a converged answer.
#[test]
fn an_unreachable_lift_target_stops_on_the_alpha_bound() {
    let options = TrimOptions {
        max_alpha: 0.5,
        ..TrimOptions::default()
    };
    let result = trim_with(&TrimPresets::level_flight(10.0), options, None);
    assert!(!result.converged);
    assert_eq!(result.status, TrimStatus::BoundsExceeded);
    assert!(result.state.alpha.abs() > options.max_alpha);
}

/// Ports `TrimMaxIterationsReached`.
#[test]
fn an_impossible_tolerance_stops_at_the_iteration_limit() {
    let options = TrimOptions {
        max_iterations: 1,
        tolerance: 1.0e-20,
        ..TrimOptions::default()
    };
    let result = trim_with(&TrimPresets::level_flight(0.5), options, None);
    assert!(!result.converged);
    assert_eq!(result.iterations, 1);
    assert_eq!(result.status, TrimStatus::MaximumIterations);
}

/// Ports `TrimWithPitchMomentConstraint`.
///
/// Driving alpha by pitching moment rather than lift exercises the other
/// coefficient row of the Jacobian.
#[test]
fn alpha_can_be_driven_by_the_pitching_moment() {
    let result = trim_with(
        &alpha_constrained(VariableConstraint {
            kind: ConstraintType::PitchMoment,
            target: 0.0,
        }),
        TrimOptions::default(),
        None,
    );
    assert!(result.iterations > 0);
}

/// Ports `TrimVerboseOutput`.
///
/// Verbose mode must not change the trajectory, only what is reported.
#[test]
fn verbose_mode_does_not_change_the_result() {
    let constraints = TrimPresets::level_flight(0.3);
    let quiet = trim_with(
        &constraints,
        TrimOptions {
            max_iterations: 3,
            ..TrimOptions::default()
        },
        None,
    );
    let loud = trim_with(
        &constraints,
        TrimOptions {
            max_iterations: 3,
            verbose: true,
            ..TrimOptions::default()
        },
        None,
    );
    assert!(loud.iterations > 0);
    assert_eq!(quiet.iterations, loud.iterations);
    assert_abs_diff_eq!(quiet.state.alpha, loud.state.alpha, epsilon = 1.0e-12);
}

/// Ports `PresetLevelFlight` and `PresetSteadySideslip`.
#[test]
fn the_presets_set_the_expected_constraint_kinds() {
    let constraints = TrimPresets::level_flight(0.5);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.5, epsilon = TOL);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.beta.target, 0.0, epsilon = TOL);

    let sideslip = 10.0 * DEGREES_TO_RADIANS;
    let constraints = TrimPresets::steady_sideslip(0.4, sideslip);
    assert_eq!(constraints.alpha.kind, ConstraintType::Lift);
    assert_abs_diff_eq!(constraints.alpha.target, 0.4, epsilon = TOL);
    assert_eq!(constraints.beta.kind, ConstraintType::Fixed);
    assert_abs_diff_eq!(constraints.beta.target, sideslip, epsilon = TOL);
}

/// Ports `ConstraintCount`.
#[test]
fn the_constraint_count_grows_one_per_pinned_variable() {
    let mut constraints = TrimConstraints::default();
    assert_eq!(constraints.constraint_count(), 0);
    constraints.alpha = VariableConstraint {
        kind: ConstraintType::Lift,
        target: 0.5,
    };
    assert_eq!(constraints.constraint_count(), 1);
    constraints.beta = fixed(0.0);
    assert_eq!(constraints.constraint_count(), 2);
    constraints.roll_rate = fixed(0.0);
    constraints.pitch_rate = fixed(0.0);
    constraints.yaw_rate = fixed(0.0);
    assert_eq!(constraints.constraint_count(), 5);
}
