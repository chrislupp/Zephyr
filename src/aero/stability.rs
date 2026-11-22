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

//! Stability derivative extraction from cached unit responses.

use crate::core::{EPSILON, Scalar, Vec3};
use crate::error::{Result, ZephyrError};
use crate::geometry::Aircraft;
use crate::solver::{FlightCondition, SolutionState, VortexLattice, VortexLatticeSolver};

use super::{
    AerodynamicResults, BodyAxisDerivatives, Coefficients, ControlDerivatives, ForceCalculator,
    StabilityDerivatives,
};

/// Computes angle and nondimensional rate derivatives.
#[derive(Clone, Copy, Debug, Default)]
pub struct StabilityDerivativeCalculator;

impl StabilityDerivativeCalculator {
    fn unit_forces(
        aircraft: &Aircraft,
        mesh: &VortexLattice,
        circulation: crate::core::VecX,
        condition: &FlightCondition,
    ) -> AerodynamicResults {
        let solution = SolutionState {
            circulation,
            is_valid: true,
            ..SolutionState::default()
        };
        let mut calculator = ForceCalculator::default();
        calculator.set_include_induced_velocity(false);
        calculator.compute_for_aircraft(aircraft, mesh, &solution, condition)
    }

    pub fn compute(
        &self,
        aircraft: &Aircraft,
        solver: &mut VortexLatticeSolver,
        condition: &FlightCondition,
    ) -> Result<StabilityDerivatives> {
        ensure_solved(solver, condition)?;
        let solution = solver.solution();
        let mesh = solver.mesh();
        let unit = [
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_u.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_v.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_w.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_p.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_q.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_r.clone(),
                condition,
            ),
        ];
        if unit.iter().any(|result| !result.is_valid) {
            return Err(ZephyrError::InvalidSolution);
        }
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        let du_da = -sin_alpha * cos_beta;
        let dw_da = cos_alpha * cos_beta;
        let du_db = -cos_alpha * sin_beta;
        let dv_db = cos_beta;
        let dw_db = -sin_alpha * sin_beta;
        let by_alpha = Coefficients::linear_combination(
            unit[0].coefficients,
            du_da,
            unit[1].coefficients,
            0.0,
            unit[2].coefficients,
            dw_da,
        );
        let by_beta = Coefficients::linear_combination(
            unit[0].coefficients,
            du_db,
            unit[1].coefficients,
            dv_db,
            unit[2].coefficients,
            dw_db,
        );
        let mut result = StabilityDerivatives {
            is_valid: true,
            reference_area: aircraft.reference_area(),
            reference_chord: aircraft.reference_chord(),
            reference_span: aircraft.reference_span(),
            velocity: 1.0,
            alpha: condition.alpha,
            beta: condition.beta,
            by_alpha,
            by_beta,
            by_roll_rate: unit[3].coefficients,
            by_pitch_rate: unit[4].coefficients,
            by_yaw_rate: unit[5].coefficients,
            ..StabilityDerivatives::default()
        };
        if result.by_alpha.lift.abs() > EPSILON {
            result.neutral_point_x = aircraft.moment_reference_point().x
                - result.by_alpha.pitch / result.by_alpha.lift * result.reference_chord;
            result.static_margin = (result.neutral_point_x - aircraft.moment_reference_point().x)
                / result.reference_chord;
        }
        for name in solver.mesh().control_surface_names() {
            let circulation = solver.control_circulation_derivative(condition, name)?;
            let response = Self::unit_forces(aircraft, mesh, circulation, condition);
            let deflection = condition
                .control_deflections
                .iter()
                .filter(|(candidate, _)| candidate == name)
                .map(|(_, value)| *value)
                .sum();
            result.control_derivatives.push(ControlDerivatives {
                name: name.clone(),
                deflection,
                gain: 1.0,
                coefficients: response.coefficients,
                body_force: response.body_force,
                hinge_moment: 0.0,
            });
        }
        Ok(result)
    }

    pub fn compute_body_axis(
        &self,
        aircraft: &Aircraft,
        solver: &mut VortexLatticeSolver,
        condition: &FlightCondition,
    ) -> Result<BodyAxisDerivatives> {
        ensure_solved(solver, condition)?;
        let solution = solver.solution();
        let mesh = solver.mesh();
        let unit = [
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_u.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_v.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_w.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_p.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_q.clone(),
                condition,
            ),
            Self::unit_forces(
                aircraft,
                mesh,
                solution.circulation_per_r.clone(),
                condition,
            ),
        ];
        let (sin_alpha, cos_alpha) = condition.alpha.sin_cos();
        let (sin_beta, cos_beta) = condition.beta.sin_cos();
        let da = [-sin_alpha * cos_beta, 0.0, cos_alpha * cos_beta];
        let db = [-cos_alpha * sin_beta, cos_beta, -sin_alpha * sin_beta];
        Ok(BodyAxisDerivatives {
            is_valid: true,
            reference_area: aircraft.reference_area(),
            reference_chord: aircraft.reference_chord(),
            reference_span: aircraft.reference_span(),
            velocity: 1.0,
            alpha: condition.alpha,
            beta: condition.beta,
            by_alpha_force: combine_vec(&unit, da, false),
            by_alpha_moment: combine_vec(&unit, da, true),
            by_beta_force: combine_vec(&unit, db, false),
            by_beta_moment: combine_vec(&unit, db, true),
            by_roll_rate_force: unit[3].body_force,
            by_roll_rate_moment: unit[3].body_moment,
            by_pitch_rate_force: unit[4].body_force,
            by_pitch_rate_moment: unit[4].body_moment,
            by_yaw_rate_force: unit[5].body_force,
            by_yaw_rate_moment: unit[5].body_moment,
        })
    }
}

/// Solve `condition` unless `solver` already holds a complete solution of it.
///
/// Solving the same condition again reproduces the stored solution bit for
/// bit, so a caller that has just solved — the trim loop does, every
/// iteration — need not pay for the right-hand side, the six unit responses
/// and the induced-velocity product a second time. A solution is reused only
/// when the unit responses it was built from are still current (moving the
/// reference point or rebuilding the AIC invalidates them without discarding
/// the solution) and it carries the per-unit circulations read below, which
/// the selective solve path does not produce.
fn ensure_solved(solver: &mut VortexLatticeSolver, condition: &FlightCondition) -> Result<()> {
    let current = solver.are_unit_responses_valid()
        && solver.check_solved_condition(condition).is_ok()
        && solver.solution().circulation_per_u.len() == solver.mesh().vortex_count();
    if !current {
        if !solver.are_unit_responses_valid() {
            solver.solve_unit_responses()?;
        }
        solver.solve(condition)?;
    }
    Ok(())
}

fn combine_vec(unit: &[AerodynamicResults; 6], weights: [Scalar; 3], moment: bool) -> Vec3 {
    let value = |index: usize| {
        if moment {
            unit[index].body_moment
        } else {
            unit[index].body_force
        }
    };
    weights[0] * value(0) + weights[1] * value(1) + weights[2] * value(2)
}
