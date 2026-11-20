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

//! Ported from the frozen v0.2.0 `adjoint_solver_test.cpp`.
//!
//! The transposed AIC solve and the circulation gradients it is seeded with.
//!
//! One case, `AdjointMatchesFiniteDiffForCL`, was `DISABLED_` in the C++ and so
//! never ran there. It is ported and enabled here: the Rust total derivative
//! includes the body-to-stability rotation term the v0.2.0 alpha expression
//! omitted, so the check the oracle could not pass now holds against a full
//! re-solve.

use approx::assert_abs_diff_eq;
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{DEGREES_TO_RADIANS, MatX, Vec3, VecX};
use zephyr_vlm::geometry::{Aircraft, Section, Surface};
use zephyr_vlm::gradients::{
    AdjointSolver, AeroOutput, TotalDerivative, force_circulation_gradient,
};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

/// The frozen fixture: a 4x6 flat wing spanning `y` in `[0, 3]`.
fn aircraft() -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(6);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 3.0, 0.0), 1.0, 0.0));

    let mut aircraft = Aircraft::new("Test");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(3.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(6.0);
    aircraft
}

fn solved(alpha_degrees: f64) -> (Aircraft, VortexLatticeSolver, FlightCondition) {
    let aircraft = aircraft();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let condition = FlightCondition {
        alpha: alpha_degrees * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    solver.solve(&condition).unwrap();
    (aircraft, solver, condition)
}

/// A deterministic pseudo-random vector, standing in for the C++ `VecX::Random`
/// so the case is reproducible rather than differently random each run.
fn pseudo_random(n: usize, seed: u64) -> VecX {
    let mut state = seed | 1;
    VecX::from_fn(n, |_, _| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ((state >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    })
}

/// Ports `ConstructFromSolverBeforeSolve`.
///
/// Without a factorized AIC there is nothing to transpose, so construction is
/// refused rather than producing an adjoint against a stale or absent matrix.
#[test]
fn an_adjoint_cannot_be_built_before_a_solve() {
    let aircraft = aircraft();
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    assert!(AdjointSolver::new(&solver).is_err());
}

/// Ports `ConstructFromSolverAfterSolve`.
#[test]
fn an_adjoint_builds_after_a_solve_with_one_unknown_per_panel() {
    let (_, solver, _) = solved(5.0);
    let adjoint = AdjointSolver::new(&solver).unwrap();
    assert!(adjoint.unknown_count() > 0);
    assert_eq!(adjoint.unknown_count(), solver.mesh().vortex_count());
}

/// Ports `SolveAdjointReturnsCorrectSize`.
#[test]
fn the_adjoint_solve_returns_one_multiplier_per_panel() {
    let (_, solver, _) = solved(0.0);
    let adjoint = AdjointSolver::new(&solver).unwrap();
    let n = adjoint.unknown_count();
    let lambda = adjoint.solve(&VecX::from_element(n, 1.0)).unwrap();
    assert_eq!(lambda.len(), n);
}

/// Ports `AdjointSolveIsCorrect`.
///
/// The defining property: `A^T lambda = rhs` to machine precision. This is what
/// makes the adjoint an exact reuse of the primal factorization rather than an
/// approximation.
#[test]
fn the_adjoint_solve_satisfies_the_transposed_system() {
    let (_, solver, _) = solved(0.0);
    let adjoint = AdjointSolver::new(&solver).unwrap();
    let n = adjoint.unknown_count();
    let rhs = pseudo_random(n, 0x2545_F491_4F6C_DD1D);
    let lambda = adjoint.solve(&rhs).unwrap();
    let residual = solver.solution_aic().transpose() * &lambda - &rhs;
    assert!(residual.norm() < 1.0e-10 * rhs.norm());
}

/// Ports `SolveMultipleAdjointProblems`.
///
/// A batched solve must give the same answer column by column as separate
/// solves would, which is what lets K outputs share one factorization.
#[test]
fn a_batched_adjoint_solve_satisfies_every_column() {
    let (_, solver, _) = solved(0.0);
    let adjoint = AdjointSolver::new(&solver).unwrap();
    let n = adjoint.unknown_count();

    let mut rhs = MatX::zeros(n, 3);
    for column in 0..3 {
        let vector = pseudo_random(n, 0x9E37_79B9_7F4A_7C15 + column as u64);
        for row in 0..n {
            rhs[(row, column)] = vector[row];
        }
    }
    let lambda = adjoint.solve_batched(&rhs).unwrap();
    assert_eq!(lambda.shape(), (n, 3));

    let transposed = solver.solution_aic().transpose();
    for column in 0..3 {
        let residual = &transposed * lambda.column(column) - rhs.column(column);
        assert!(residual.norm() < 1.0e-10 * rhs.column(column).norm());
    }
}

/// Ports `ForceGradientsHaveCorrectSize` and `ForceGradientsAreNonZero`.
#[test]
fn every_force_circulation_gradient_is_sized_and_nonzero() {
    let (_, solver, condition) = solved(5.0);
    let n = solver.solution().circulation.len();
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
    ] {
        let gradient = force_circulation_gradient(&solver, &condition, output);
        assert_eq!(gradient.len(), n);
    }
    // Lift and pitching moment must actually respond to circulation; the
    // lateral ones are legitimately zero on a symmetric wing.
    assert!(force_circulation_gradient(&solver, &condition, AeroOutput::Lift).norm() > 0.0);
    assert!(force_circulation_gradient(&solver, &condition, AeroOutput::Pitch).norm() > 0.0);
}

/// Ports `CLGradientIsPositive`.
///
/// More circulation anywhere on a lifting wing means more lift, so the
/// overwhelming majority of entries are positive.
#[test]
fn most_lift_circulation_sensitivities_are_positive() {
    let (_, solver, condition) = solved(5.0);
    let gradient = force_circulation_gradient(&solver, &condition, AeroOutput::Lift);
    let positive = gradient.iter().filter(|value| **value > 0.0).count();
    let fraction = positive as f64 / gradient.len() as f64;
    assert!(fraction > 0.8, "only {fraction} of entries were positive");
}

/// Ports `AdjointIdentityTest`.
#[test]
fn the_lift_adjoint_multiplier_is_nontrivial() {
    let (_, solver, condition) = solved(5.0);
    let gradient = force_circulation_gradient(&solver, &condition, AeroOutput::Lift);
    let lambda = AdjointSolver::new(&solver)
        .unwrap()
        .solve(&gradient)
        .unwrap();
    assert_eq!(lambda.len(), solver.solution().circulation.len());
    assert!(lambda.norm() > 1.0e-10);
}

/// Ports `AdjointValidationTest.AdjointMatchesFiniteDiffForCL`, which was
/// `DISABLED_` upstream and never ran.
///
/// The total derivative is the explicit term plus the implicit one carried by
/// the adjoint. v0.2.0's explicit alpha expression omitted the
/// body-to-stability rotation, which is presumably why the case was disabled;
/// the Rust includes it, so the adjoint now agrees with a full re-solve far
/// inside the 10% band the frozen case allowed.
#[test]
fn the_adjoint_lift_derivative_matches_a_full_resolve() {
    // The smaller `AdjointValidationTest` fixture.
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(3);
    wing.set_spanwise_panel_count(4);
    wing.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Validation");
    aircraft.add_surface(wing);
    aircraft.set_reference_area(2.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);

    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let adjoint = TotalDerivative::compute(&solver, &condition, AeroOutput::Lift).unwrap();

    let step = 1.0e-6;
    let lift_at = |alpha: f64| -> f64 {
        // `FlightCondition` owns its control-deflection list, so the base is
        // cloned rather than moved out of the enclosing scope.
        let perturbed = FlightCondition {
            alpha,
            ..condition.clone()
        };
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft).unwrap();
        let solution = solver.solve(&perturbed).unwrap();
        ForceCalculator::new()
            .compute_for_aircraft(&aircraft, solver.mesh(), &solution, &perturbed)
            .coefficients
            .lift
    };
    let difference =
        (lift_at(condition.alpha + step) - lift_at(condition.alpha - step)) / (2.0 * step);

    // The frozen case allowed 10% plus an absolute 0.01; the analytical
    // derivative lands within the difference's own truncation error.
    assert_abs_diff_eq!(
        adjoint.by_alpha,
        difference,
        epsilon = 1.0e-6 * difference.abs()
    );
}
