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

//! Damped Newton trim of a prescribed, possibly deformed, direct mesh.

use crate::core::{MatX, VecX};
use crate::error::{Result, ZephyrError};
use crate::geometry::MeshInput;
use crate::solver::{FlightCondition, VortexLatticeSolver};

use super::evaluation::linearize_trim_states;
use super::{MeshTrimProblem, aerodynamic_values};

/// Numerical diagnostics, with residuals in their physical quantity units.
#[derive(Clone, Debug)]
pub struct TrimDiagnostics {
    /// Accepted Newton steps.
    pub iterations: usize,
    /// Final unscaled residuals, in balance-row order.
    pub residuals: Vec<f64>,
    /// Singular-value ratio of the scaled final Jacobian.
    pub singular_value_ratio: f64,
}

/// Converged mesh trim and the linearization used for exact implicit sensitivities.
#[derive(Clone, Debug)]
pub struct MeshTrimSolution {
    /// Solved operating condition, including realized surface commands.
    pub condition: FlightCondition,
    /// Matching converged aerodynamic state and cached influence factorization.
    pub solver: VortexLatticeSolver,
    /// Solved independent states.
    pub states: VecX,
    /// Coefficient then dimensional resultant values.
    pub values: [f64; 12],
    /// Configuration of the converged residual problem.
    pub problem: MeshTrimProblem,
    /// Output-state derivatives selected by the balance equations. Other rows
    /// may be zero; `compute_partials` selects the requested functions as well.
    pub state_partials: MatX,
    /// Scaled residual Jacobian with respect to physical state values.
    pub jacobian: MatX,
    /// Convergence and conditioning diagnostics.
    pub diagnostics: TrimDiagnostics,
    pub(crate) mesh: MeshInput,
    pub(crate) factorization: nalgebra::linalg::LU<f64, nalgebra::Dyn, nalgebra::Dyn>,
}

/// Internal trim root finder, independent of MDO transport.
#[derive(Clone, Copy, Debug, Default)]
pub struct MeshTrimSolver;

fn failure(message: &str, z: &VecX, residual: &VecX, problem: &MeshTrimProblem) -> ZephyrError {
    let rows: Vec<_> = problem
        .balances
        .iter()
        .zip(residual.iter())
        .map(|(b, r)| (b.quantity.name(), r * b.scale))
        .collect();
    ZephyrError::TrimFailed(format!(
        "{message}; states={:?}; physical residuals={rows:?}",
        z.as_slice()
    ))
}

impl MeshTrimSolver {
    /// Converge a square trim system for prescribed mesh, targets, and external loads.
    ///
    /// Initial guesses are configuration. No derivatives of solver iterations or
    /// guesses enter the converged explicit mapping.
    pub fn solve(
        &self,
        mesh: &MeshInput,
        base: &FlightCondition,
        problem: &MeshTrimProblem,
        targets: &[f64],
        external: &[f64; 6],
    ) -> Result<MeshTrimSolution> {
        problem.validate(mesh)?;
        problem.residual(&[0.0; 12], targets, external)?;
        let mut solver = VortexLatticeSolver::new();
        solver.set_mesh(mesh)?;
        let mut z = VecX::from_iterator(
            problem.variables.len(),
            problem.variables.iter().map(|v| v.initial_guess),
        );
        let mut condition = problem.condition(base, &z)?;
        solver.solve(&condition)?;
        let mut values = aerodynamic_values(mesh, &solver, &condition);
        let mut residual = problem.residual(&values, targets, external)?;
        let initial = residual.map(f64::abs);
        let quantities: Vec<_> = problem.balances.iter().map(|b| b.quantity).collect();
        for iteration in 0..=problem.options.max_iterations {
            let aero = linearize_trim_states(
                mesh,
                &solver,
                &condition,
                &quantities,
                problem.variables.iter().any(|v| v.name == "beta"),
            )?;
            let state_partials = problem.state_partials(&aero);
            let jacobian = problem.residual_jacobian(&state_partials);
            let scaled = MatX::from_fn(jacobian.nrows(), jacobian.ncols(), |r, c| {
                jacobian[(r, c)] * problem.variables[c].scale
            });
            if scaled.iter().any(|v| !v.is_finite()) || residual.iter().any(|v| !v.is_finite()) {
                return Err(failure(
                    "nonfinite trim linearization",
                    &z,
                    &residual,
                    problem,
                ));
            }
            let sv = scaled.clone().svd(false, false).singular_values;
            let ratio = if sv.max() > 0.0 {
                sv.min() / sv.max()
            } else {
                0.0
            };
            if ratio <= problem.options.rank_tolerance {
                return Err(failure(
                    "singular or ineffective trim variables",
                    &z,
                    &residual,
                    problem,
                ));
            }
            if residual.iter().zip(initial.iter()).all(|(r, first)| {
                r.abs()
                    <= problem.options.absolute_tolerance
                        + problem.options.relative_tolerance * first
            }) {
                let diagnostics = TrimDiagnostics {
                    iterations: iteration,
                    residuals: residual
                        .iter()
                        .zip(&problem.balances)
                        .map(|(r, b)| r * b.scale)
                        .collect(),
                    singular_value_ratio: ratio,
                };
                return Ok(MeshTrimSolution {
                    condition,
                    solver,
                    states: z,
                    values,
                    problem: problem.clone(),
                    state_partials,
                    mesh: mesh.clone(),
                    factorization: jacobian.clone().lu(),
                    jacobian,
                    diagnostics,
                });
            }
            if iteration == problem.options.max_iterations {
                break;
            }
            let scaled_step = scaled
                .lu()
                .solve(&(-&residual))
                .ok_or_else(|| failure("singular Newton system", &z, &residual, problem))?;
            let step = VecX::from_fn(z.len(), |i, _| scaled_step[i] * problem.variables[i].scale);
            let norm = residual.norm();
            let mut accepted = None;
            let mut fraction = 1.0;
            for _ in 0..problem.options.max_backtracks {
                let trial_z = &z + fraction * &step;
                if let Ok(trial_condition) = problem.condition(base, &trial_z)
                    && solver.solve(&trial_condition).is_ok()
                {
                    let trial_values = aerodynamic_values(mesh, &solver, &trial_condition);
                    let trial_residual = problem.residual(&trial_values, targets, external)?;
                    if trial_residual.iter().all(|v| v.is_finite())
                        && trial_residual.norm() < (1.0 - 1e-4 * fraction) * norm
                    {
                        accepted = Some((trial_z, trial_condition, trial_values, trial_residual));
                        break;
                    }
                }
                fraction *= 0.5;
            }
            let Some((next_z, next_condition, next_values, next_residual)) = accepted else {
                return Err(failure(
                    "line search could not find an admissible improving step",
                    &z,
                    &residual,
                    problem,
                ));
            };
            z = next_z;
            condition = next_condition;
            values = next_values;
            residual = next_residual;
        }
        Err(failure("iteration limit exceeded", &z, &residual, problem))
    }
}
