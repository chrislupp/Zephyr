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

//! Newton–Raphson trim iteration.

use crate::aero::{
    AerodynamicResults, Coefficients, ForceCalculator, StabilityDerivativeCalculator,
    StabilityDerivatives,
};
use crate::core::{MatX, Scalar, Vec3, VecX};
use crate::error::Result;
use crate::geometry::Aircraft;
use crate::gradients::{
    AdjointSolver, AeroOutput, force_alpha_gradient, force_beta_gradient,
    force_circulation_gradient,
};
use crate::solver::{FlightCondition, VortexLatticeSolver};

use super::{ConstraintType, OperatingVariables, TrimConstraints, VariableConstraint};

/// Terminal state of a trim attempt.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum TrimStatus {
    #[default]
    NotStarted,
    Converged,
    NoConstraints,
    AerodynamicFailure,
    SingularJacobian,
    BoundsExceeded,
    MaximumIterations,
}

/// Solver limits and convergence controls.
#[derive(Clone, Copy, Debug)]
pub struct TrimOptions {
    pub max_iterations: usize,
    pub tolerance: Scalar,
    pub max_alpha: Scalar,
    pub max_beta: Scalar,
    pub max_rate: Scalar,
    pub relaxation: Scalar,
    pub verbose: bool,
}

impl Default for TrimOptions {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            tolerance: 1.0e-6,
            max_alpha: 1.5,
            max_beta: 1.5,
            max_rate: 5.0,
            relaxation: 1.0,
            verbose: false,
        }
    }
}

/// Normal, non-exceptional result of a trim attempt.
#[derive(Clone, Debug, Default)]
pub struct TrimResult {
    pub status: TrimStatus,
    pub converged: bool,
    pub iterations: usize,
    pub residual_norm: Scalar,
    pub state: OperatingVariables,
    pub aerodynamic_results: AerodynamicResults,
    pub derivatives: StabilityDerivatives,
    pub diagnostic: String,
}

/// Newton trim driver.
#[derive(Clone, Debug, Default)]
pub struct TrimSolver {
    options: TrimOptions,
}

impl TrimSolver {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_options(&mut self, options: TrimOptions) {
        self.options = options;
    }

    #[must_use]
    pub fn options(&self) -> TrimOptions {
        self.options
    }

    #[must_use]
    pub fn solve(
        &self,
        aircraft: &Aircraft,
        solver: &mut VortexLatticeSolver,
        constraints: &TrimConstraints,
        initial: Option<&OperatingVariables>,
    ) -> TrimResult {
        let mut result = TrimResult {
            state: initial.cloned().unwrap_or_default(),
            ..TrimResult::default()
        };
        let count = constraints.constraint_count();
        if count == 0 {
            result.status = TrimStatus::NoConstraints;
            result.diagnostic = "no constraints specified".into();
            return result;
        }
        if solver.mesh().vortex_count() == 0
            && let Err(error) = solver.set_aircraft(aircraft)
        {
            return aerodynamic_failure(result, error.to_string());
        }
        let force_calculator = ForceCalculator::default();
        let derivative_calculator = StabilityDerivativeCalculator;
        for iteration in 0..self.options.max_iterations {
            result.iterations = iteration + 1;
            let condition = FlightCondition {
                alpha: result.state.alpha,
                beta: result.state.beta,
                roll_rate: result.state.roll_rate,
                pitch_rate: result.state.pitch_rate,
                yaw_rate: result.state.yaw_rate,
                control_deflections: result.state.control_deflections.clone(),
                ..FlightCondition::default()
            };
            let solution = match solver.solve(&condition) {
                Ok(solution) => solution,
                Err(error) => {
                    return aerodynamic_failure(result, error.to_string());
                }
            };
            result.aerodynamic_results = force_calculator.compute_for_aircraft(
                aircraft,
                solver.mesh(),
                &solution,
                &condition,
            );
            if !result.aerodynamic_results.is_valid {
                return aerodynamic_failure(result, "force computation failed".into());
            }
            // The reported derivatives keep their frozen-induced-velocity
            // convention; `compute` reuses the solution just obtained rather
            // than solving the same condition again.
            result.derivatives = match derivative_calculator.compute(aircraft, solver, &condition) {
                Ok(derivatives) => derivatives,
                Err(error) => {
                    return aerodynamic_failure(result, error.to_string());
                }
            };
            let residual = trim_residual(&result.aerodynamic_results, &result.state, constraints);
            result.residual_norm = residual.norm();
            if self.options.verbose {
                eprintln!(
                    "trim iteration {}: residual={:.6e}, alpha={:.6} rad",
                    result.iterations, result.residual_norm, result.state.alpha
                );
            }
            if result.residual_norm < self.options.tolerance {
                result.status = TrimStatus::Converged;
                result.converged = true;
                return result;
            }
            // The Newton step uses exact total derivatives of the residual,
            // not the reported ones: those hold the induced velocity frozen
            // and, for alpha and beta, also leave out the turning of the
            // freestream and of the stability axes, so a step built on them
            // converges only linearly.
            let sensitivities = match TrimSensitivities::compute(solver, &condition, constraints) {
                Ok(sensitivities) => sensitivities,
                Err(error) => {
                    return aerodynamic_failure(result, error.to_string());
                }
            };
            let jacobian = trim_jacobian(&sensitivities, constraints);
            let Some(mut delta) = jacobian.lu().solve(&(-residual)) else {
                result.status = TrimStatus::SingularJacobian;
                result.diagnostic = "trim Jacobian is singular".into();
                return result;
            };
            delta *= self.options.relaxation;
            update_state(&mut result.state, &delta, constraints);
            if !self.in_bounds(&result.state) {
                result.status = TrimStatus::BoundsExceeded;
                result.diagnostic = "operating variables exceeded configured bounds".into();
                return result;
            }
        }
        result.status = TrimStatus::MaximumIterations;
        result.diagnostic = "maximum iterations reached".into();
        result
    }

    fn in_bounds(&self, state: &OperatingVariables) -> bool {
        state.alpha.abs() <= self.options.max_alpha
            && state.beta.abs() <= self.options.max_beta
            && state.roll_rate.abs() <= self.options.max_rate
            && state.pitch_rate.abs() <= self.options.max_rate
            && state.yaw_rate.abs() <= self.options.max_rate
    }
}

fn aerodynamic_failure(mut result: TrimResult, diagnostic: String) -> TrimResult {
    result.status = TrimStatus::AerodynamicFailure;
    result.diagnostic = diagnostic;
    result
}

fn active(constraint: VariableConstraint) -> bool {
    constraint.kind != ConstraintType::Free
}

/// Coefficient an equation of the given kind drives, if any.
fn output_of(kind: ConstraintType) -> Option<AeroOutput> {
    match kind {
        ConstraintType::Lift => Some(AeroOutput::Lift),
        ConstraintType::Side => Some(AeroOutput::Side),
        ConstraintType::RollMoment => Some(AeroOutput::Roll),
        ConstraintType::PitchMoment => Some(AeroOutput::Pitch),
        ConstraintType::YawMoment => Some(AeroOutput::Yaw),
        ConstraintType::Free | ConstraintType::Fixed => None,
    }
}

fn select(coefficients: Coefficients, output: AeroOutput) -> Scalar {
    match output {
        AeroOutput::Lift => coefficients.lift,
        AeroOutput::Drag => coefficients.drag,
        AeroOutput::Side => coefficients.side,
        AeroOutput::Roll => coefficients.roll,
        AeroOutput::Pitch => coefficients.pitch,
        AeroOutput::Yaw => coefficients.yaw,
        AeroOutput::InducedDrag => 0.0,
    }
}

fn select_mut(coefficients: &mut Coefficients, output: AeroOutput) -> Option<&mut Scalar> {
    match output {
        AeroOutput::Lift => Some(&mut coefficients.lift),
        AeroOutput::Drag => Some(&mut coefficients.drag),
        AeroOutput::Side => Some(&mut coefficients.side),
        AeroOutput::Roll => Some(&mut coefficients.roll),
        AeroOutput::Pitch => Some(&mut coefficients.pitch),
        AeroOutput::Yaw => Some(&mut coefficients.yaw),
        AeroOutput::InducedDrag => None,
    }
}

/// Index of a coefficient in the six-component explicit-gradient arrays.
fn explicit_index(output: AeroOutput) -> usize {
    match output {
        AeroOutput::Lift | AeroOutput::InducedDrag => 0,
        AeroOutput::Drag => 1,
        AeroOutput::Side => 2,
        AeroOutput::Roll => 3,
        AeroOutput::Pitch => 4,
        AeroOutput::Yaw => 5,
    }
}

/// Exact total derivatives of the trimmed coefficients at a solved condition.
///
/// Only the coefficients some equation drives and the variables some
/// constraint activates are computed; every other entry stays zero. The
/// entries are the same adjoint contractions `GradientCalculator` uses, which
/// differentiate the reported near-field forces including the induced
/// velocity, so a Newton step on them converges quadratically.
#[derive(Clone, Debug, Default)]
struct TrimSensitivities {
    /// Derivatives by alpha, beta, roll, pitch, and yaw rate.
    states: [Coefficients; 5],
    /// Derivatives by each constrained control; `None` for a name the lattice
    /// does not declare.
    controls: Vec<Option<Coefficients>>,
}

impl TrimSensitivities {
    fn compute(
        solver: &VortexLatticeSolver,
        condition: &FlightCondition,
        constraints: &TrimConstraints,
    ) -> Result<Self> {
        let states = state_constraints(constraints);
        let known = |name: &str| {
            solver
                .mesh()
                .control_surface_names()
                .iter()
                .any(|candidate| candidate == name)
        };
        let mut result = Self {
            controls: constraints
                .controls
                .iter()
                .map(|(name, _)| known(name).then(Coefficients::default))
                .collect(),
            ..Self::default()
        };
        let mut outputs: Vec<AeroOutput> = Vec::new();
        for kind in states
            .iter()
            .chain(
                constraints
                    .controls
                    .iter()
                    .map(|(_, constraint)| constraint),
            )
            .filter(|constraint| active(**constraint))
            .filter_map(|constraint| output_of(constraint.kind))
        {
            if !outputs.contains(&kind) {
                outputs.push(kind);
            }
        }
        if outputs.is_empty() {
            return Ok(result);
        }
        let mesh = solver.mesh();
        let solution = solver.solution();
        let reference = solver.reference_quantities();
        let n = mesh.vortex_count();
        let mut seeds = MatX::zeros(n, outputs.len());
        for (column, output) in outputs.iter().enumerate() {
            seeds.set_column(
                column,
                &force_circulation_gradient(solver, condition, *output),
            );
        }
        let adjoints = AdjointSolver::new(solver)?.solve_batched(&seeds)?;

        // Freestream sensitivities of the boundary-condition residual.
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        let onset = |velocity: Vec3| {
            VecX::from_fn(n, |row, _| {
                if mesh.vortices()[row].affected_by_freestream {
                    -solver.linearized_normals()[row].dot(&velocity)
                } else {
                    0.0
                }
            })
        };
        let explicit_alpha = active(states[0]).then(|| {
            (
                force_alpha_gradient(mesh, solution, condition, reference),
                onset(Vec3::new(-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta)),
            )
        });
        let explicit_beta = active(states[1]).then(|| {
            (
                force_beta_gradient(mesh, solution, condition, reference),
                onset(Vec3::new(
                    -cos_alpha * sin_beta,
                    cos_beta,
                    -sin_alpha * sin_beta,
                )),
            )
        });
        let rate_circulations = [
            &solution.circulation_per_p,
            &solution.circulation_per_q,
            &solution.circulation_per_r,
        ];
        let control_residuals = constraints
            .controls
            .iter()
            .zip(&result.controls)
            .map(|((name, constraint), known)| {
                if active(*constraint) && known.is_some() {
                    solver
                        .control_residual_derivative(condition, name)
                        .map(Some)
                } else {
                    Ok(None)
                }
            })
            .collect::<Result<Vec<_>>>()?;

        for (column, output) in outputs.iter().copied().enumerate() {
            let seed = seeds.column(column);
            let adjoint = adjoints.column(column);
            let index = explicit_index(output);
            let mut values = [0.0; 5];
            if let Some((explicit, rhs)) = &explicit_alpha {
                values[0] = explicit[index] + adjoint.dot(rhs);
            }
            if let Some((explicit, rhs)) = &explicit_beta {
                values[1] = explicit[index] + adjoint.dot(rhs);
            }
            for (slot, circulation) in rate_circulations.iter().enumerate() {
                if active(states[slot + 2]) {
                    values[slot + 2] = seed.dot(*circulation);
                }
            }
            for (state, value) in result.states.iter_mut().zip(values) {
                if let Some(entry) = select_mut(state, output) {
                    *entry = value;
                }
            }
            for (control, rhs) in result.controls.iter_mut().zip(&control_residuals) {
                if let (Some(control), Some(rhs)) = (control.as_mut(), rhs)
                    && let Some(entry) = select_mut(control, output)
                {
                    *entry = adjoint.dot(rhs);
                }
            }
        }
        Ok(result)
    }
}

fn state_constraints(constraints: &TrimConstraints) -> [VariableConstraint; 5] {
    [
        constraints.alpha,
        constraints.beta,
        constraints.roll_rate,
        constraints.pitch_rate,
        constraints.yaw_rate,
    ]
}

/// Active equations in row order, each with the coefficient it drives (`None`
/// for a fixed-value equation). Row `k` belongs to active variable `k`.
fn equations(constraints: &TrimConstraints) -> Vec<(VariableConstraint, Option<AeroOutput>)> {
    state_constraints(constraints)
        .into_iter()
        .chain(
            constraints
                .controls
                .iter()
                .map(|(_, constraint)| *constraint),
        )
        .filter(|constraint| active(*constraint))
        .map(|constraint| (constraint, output_of(constraint.kind)))
        .collect()
}

/// Trim residual `r`, one entry per active equation.
fn trim_residual(
    aero: &AerodynamicResults,
    state: &OperatingVariables,
    constraints: &TrimConstraints,
) -> VecX {
    let deflection = |name: &str| {
        state
            .control_deflections
            .iter()
            .find(|(candidate, _)| candidate == name)
            .map_or(0.0, |(_, value)| *value)
    };
    let values = [
        state.alpha,
        state.beta,
        state.roll_rate,
        state.pitch_rate,
        state.yaw_rate,
    ]
    .into_iter()
    .zip(state_constraints(constraints))
    .chain(
        constraints
            .controls
            .iter()
            .map(|(name, constraint)| (deflection(name), *constraint)),
    )
    .filter(|(_, constraint)| active(*constraint))
    .map(|(value, _)| value);
    VecX::from_iterator(
        constraints.constraint_count(),
        equations(constraints)
            .into_iter()
            .zip(values)
            .map(|((constraint, output), value)| {
                output.map_or(value, |output| select(aero.coefficients, output)) - constraint.target
            }),
    )
}

/// Newton Jacobian `dr/dx` over the active variables.
///
/// Every coefficient equation is differentiated with respect to every active
/// variable, state and control alike. Assembling each equation only against
/// variables of its own kind left out the control column of a state equation
/// and the state columns of a control equation; for a lift-and-elevator trim
/// that dropped `dCL/d(elevator)` and `dCm/d(alpha)` and turned Newton into a
/// fixed-point iteration contracting by about 0.3 per step.
fn trim_jacobian(sensitivities: &TrimSensitivities, constraints: &TrimConstraints) -> MatX {
    // Derivatives of every coefficient by each active variable, in column
    // order; `None` marks a control the lattice does not declare.
    let columns: Vec<Option<Coefficients>> = state_constraints(constraints)
        .into_iter()
        .zip(sensitivities.states.iter().copied().map(Some))
        .chain(
            constraints
                .controls
                .iter()
                .map(|(_, constraint)| *constraint)
                .zip(sensitivities.controls.iter().copied()),
        )
        .filter(|(constraint, _)| active(*constraint))
        .map(|(_, derivative)| derivative)
        .collect();
    let count = columns.len();
    let mut jacobian = MatX::zeros(count, count);
    // Equation `k` belongs to variable `k`, so its own column is the diagonal.
    for (row, (_, output)) in equations(constraints).into_iter().enumerate() {
        match output {
            Some(output) => {
                for (column, derivative) in columns.iter().enumerate() {
                    if let Some(derivative) = derivative {
                        jacobian[(row, column)] = select(*derivative, output);
                    }
                }
                if columns[row].is_none() {
                    // A control the lattice does not declare has no
                    // derivative. The unit placeholder lets the step be taken
                    // so the next aerodynamic solve reports the unknown name.
                    jacobian[(row, row)] = 1.0;
                }
            }
            None => jacobian[(row, row)] = 1.0,
        }
    }
    jacobian
}

fn update_state(state: &mut OperatingVariables, delta: &VecX, constraints: &TrimConstraints) {
    let mut index = 0;
    for (constraint, value) in [
        (constraints.alpha, &mut state.alpha),
        (constraints.beta, &mut state.beta),
        (constraints.roll_rate, &mut state.roll_rate),
        (constraints.pitch_rate, &mut state.pitch_rate),
        (constraints.yaw_rate, &mut state.yaw_rate),
    ] {
        if active(constraint) {
            *value += delta[index];
            index += 1;
        }
    }
    for (name, constraint) in &constraints.controls {
        if !active(*constraint) {
            continue;
        }
        if let Some((_, value)) = state
            .control_deflections
            .iter_mut()
            .find(|(candidate, _)| candidate == name)
        {
            *value += delta[index];
        } else {
            state.control_deflections.push((name.clone(), delta[index]));
        }
        index += 1;
    }
}
