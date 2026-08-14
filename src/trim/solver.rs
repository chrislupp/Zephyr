//! Newton–Raphson trim iteration.

use crate::aero::{
    AerodynamicResults, ForceCalculator, StabilityDerivativeCalculator, StabilityDerivatives,
};
use crate::core::{MatX, Scalar, VecX};
use crate::geometry::Aircraft;
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
            result.derivatives = match derivative_calculator.compute(aircraft, solver, &condition) {
                Ok(derivatives) => derivatives,
                Err(error) => {
                    return aerodynamic_failure(result, error.to_string());
                }
            };
            let (jacobian, residual) = build_system(
                &result.aerodynamic_results,
                &result.derivatives,
                &result.state,
                constraints,
            );
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

fn variable_index(constraints: &TrimConstraints, target: usize) -> usize {
    [
        constraints.alpha,
        constraints.beta,
        constraints.roll_rate,
        constraints.pitch_rate,
        constraints.yaw_rate,
    ][..target]
        .iter()
        .filter(|constraint| active(**constraint))
        .count()
}

fn build_system(
    aero: &AerodynamicResults,
    derivatives: &StabilityDerivatives,
    state: &OperatingVariables,
    constraints: &TrimConstraints,
) -> (MatX, VecX) {
    let count = constraints.constraint_count();
    let mut jacobian = MatX::zeros(count, count);
    let mut residual = VecX::zeros(count);
    let mut row = 0;
    let base = [
        constraints.alpha,
        constraints.beta,
        constraints.roll_rate,
        constraints.pitch_rate,
        constraints.yaw_rate,
    ];
    let state_values = [
        state.alpha,
        state.beta,
        state.roll_rate,
        state.pitch_rate,
        state.yaw_rate,
    ];
    for (slot, constraint) in base.iter().copied().enumerate() {
        if !active(constraint) {
            continue;
        }
        match constraint.kind {
            ConstraintType::Fixed => {
                residual[row] = state_values[slot] - constraint.target;
                jacobian[(row, variable_index(constraints, slot))] = 1.0;
            }
            ConstraintType::Lift => {
                residual[row] = aero.coefficients.lift - constraint.target;
                coefficient_row(
                    &mut jacobian,
                    row,
                    constraints,
                    derivatives,
                    |coefficient| coefficient.lift,
                );
            }
            ConstraintType::Side => {
                residual[row] = aero.coefficients.side - constraint.target;
                coefficient_row(
                    &mut jacobian,
                    row,
                    constraints,
                    derivatives,
                    |coefficient| coefficient.side,
                );
            }
            ConstraintType::RollMoment => {
                residual[row] = aero.coefficients.roll - constraint.target;
                coefficient_row(
                    &mut jacobian,
                    row,
                    constraints,
                    derivatives,
                    |coefficient| coefficient.roll,
                );
            }
            ConstraintType::PitchMoment => {
                residual[row] = aero.coefficients.pitch - constraint.target;
                coefficient_row(
                    &mut jacobian,
                    row,
                    constraints,
                    derivatives,
                    |coefficient| coefficient.pitch,
                );
            }
            ConstraintType::YawMoment => {
                residual[row] = aero.coefficients.yaw - constraint.target;
                coefficient_row(
                    &mut jacobian,
                    row,
                    constraints,
                    derivatives,
                    |coefficient| coefficient.yaw,
                );
            }
            ConstraintType::Free => {}
        }
        row += 1;
    }
    let control_offset = base
        .iter()
        .filter(|constraint| active(**constraint))
        .count();
    for (control_slot, (name, constraint)) in constraints.controls.iter().enumerate() {
        if !active(*constraint) {
            continue;
        }
        let column = control_offset
            + constraints.controls[..control_slot]
                .iter()
                .filter(|(_, prior)| active(*prior))
                .count();
        let deflection = state
            .control_deflections
            .iter()
            .find(|(candidate, _)| candidate == name)
            .map_or(0.0, |(_, value)| *value);
        match constraint.kind {
            ConstraintType::Fixed => {
                residual[row] = deflection - constraint.target;
                jacobian[(row, column)] = 1.0;
            }
            kind => {
                let output = match kind {
                    ConstraintType::Lift => aero.coefficients.lift,
                    ConstraintType::Side => aero.coefficients.side,
                    ConstraintType::RollMoment => aero.coefficients.roll,
                    ConstraintType::PitchMoment => aero.coefficients.pitch,
                    ConstraintType::YawMoment => aero.coefficients.yaw,
                    ConstraintType::Free | ConstraintType::Fixed => 0.0,
                };
                residual[row] = output - constraint.target;
                let derivative = derivatives
                    .control_derivatives
                    .iter()
                    .find(|control| control.name == *name)
                    .map_or(1.0, |control| match kind {
                        ConstraintType::Lift => control.coefficients.lift,
                        ConstraintType::Side => control.coefficients.side,
                        ConstraintType::RollMoment => control.coefficients.roll,
                        ConstraintType::PitchMoment => control.coefficients.pitch,
                        ConstraintType::YawMoment => control.coefficients.yaw,
                        ConstraintType::Free | ConstraintType::Fixed => 1.0,
                    });
                jacobian[(row, column)] = derivative;
            }
        }
        row += 1;
    }
    (jacobian, residual)
}

fn coefficient_row(
    jacobian: &mut MatX,
    row: usize,
    constraints: &TrimConstraints,
    derivatives: &StabilityDerivatives,
    select: impl Fn(crate::aero::Coefficients) -> Scalar,
) {
    let values = [
        select(derivatives.by_alpha),
        select(derivatives.by_beta),
        select(derivatives.by_roll_rate),
        select(derivatives.by_pitch_rate),
        select(derivatives.by_yaw_rate),
    ];
    let base = [
        constraints.alpha,
        constraints.beta,
        constraints.roll_rate,
        constraints.pitch_rate,
        constraints.yaw_rate,
    ];
    for (slot, constraint) in base.iter().enumerate() {
        if active(*constraint) {
            jacobian[(row, variable_index(constraints, slot))] = values[slot];
        }
    }
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
