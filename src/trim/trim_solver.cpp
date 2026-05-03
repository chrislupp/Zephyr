/*
    Project Pegasus
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/trim/trim_solver.h"

#include "pegasus/aero/force_calculator.h"

#include <cmath>
#include <iostream>

namespace pegasus {

namespace {

/// Variable type indices for building the system
constexpr int kVarAlpha = 0;
constexpr int kVarBeta = 1;
constexpr int kVarRollRate = 2;
constexpr int kVarPitchRate = 3;
constexpr int kVarYawRate = 4;
constexpr int kVarControlBase = 5;

}  // namespace

TrimResult TrimSolver::Solve(const Aircraft& aircraft,
                             VortexLatticeSolver& solver,
                             const TrimConstraints& constraints,
                             const OperatingVariables& initial_state) {
    TrimResult result;
    result.state = initial_state;

    // Ensure we have constraints
    Index num_constraints = constraints.GetNumConstraints();
    if (num_constraints == 0) {
        result.error_message = "No constraints specified";
        return result;
    }

    // Ensure aircraft is set up
    if (solver.GetMesh().GetNumVortices() == 0) {
        solver.SetAircraft(aircraft);
    }

    ForceCalculator force_calculator;
    StabilityDerivativeCalculator deriv_calculator;

    // Newton iteration
    for (int iter = 0; iter < options_.max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Build flight condition from current state
        FlightCondition condition;
        condition.alpha = result.state.alpha;
        condition.beta = result.state.beta;
        condition.roll_rate = result.state.roll_rate;
        condition.pitch_rate = result.state.pitch_rate;
        condition.yaw_rate = result.state.yaw_rate;
        condition.control_deflections = result.state.control_deflections;

        // Solve for circulation
        SolutionState solution = solver.Solve(condition);
        if (!solution.is_valid) {
            result.error_message = "VLM solution failed";
            return result;
        }

        // Compute aerodynamic forces
        result.aero = force_calculator.ComputeForces(aircraft, solver.GetMesh(),
                                                     solution, condition);
        if (!result.aero.is_valid) {
            result.error_message = "Force computation failed";
            return result;
        }

        // Compute stability derivatives
        result.derivatives =
            deriv_calculator.ComputeDerivatives(aircraft, solver, condition);

        // Build Newton system
        MatX jacobian;
        VecX residual;
        BuildSystem(result.aero, result.derivatives, result.state, constraints,
                    jacobian, residual);

        // Check convergence
        result.residual_norm = residual.norm();
        if (options_.verbose) {
            std::cout << "Trim iter " << iter + 1
                      << ": residual = " << result.residual_norm
                      << ", alpha = " << result.state.alpha / kDegreesToRadians
                      << " deg\n";
        }

        if (result.residual_norm < options_.tolerance) {
            result.converged = true;
            return result;
        }

        // Solve Newton system: J * delta = -residual
        Eigen::PartialPivLU<MatX> lu(jacobian);
        VecX delta = lu.solve(-residual);

        // Apply relaxation
        delta *= options_.relaxation;

        // Update state
        if (!UpdateState(result.state, delta, constraints)) {
            result.error_message = "State update exceeded bounds";
            return result;
        }

        // Check bounds
        if (!CheckBounds(result.state)) {
            result.error_message = "Operating variables out of bounds";
            return result;
        }
    }

    result.error_message = "Maximum iterations reached";
    return result;
}

void TrimSolver::BuildSystem(const AerodynamicResults& aero,
                             const StabilityDerivatives& derivs,
                             const OperatingVariables& state,
                             const TrimConstraints& constraints, MatX& jacobian,
                             VecX& residual) const {
    Index n = constraints.GetNumConstraints();
    jacobian.setZero(n, n);
    residual.setZero(n);

    Index row = 0;

    // Process alpha constraint
    if (constraints.alpha.type != ConstraintType::kFree) {
        switch (constraints.alpha.type) {
            case ConstraintType::kFixed:
                residual(row) = state.alpha - constraints.alpha.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kCL:
                residual(row) = aero.CL - constraints.alpha.target;
                // Jacobian entries for CL constraint
                jacobian(row, GetVariableIndex(constraints, kVarAlpha)) =
                    derivs.CL_alpha;
                if (constraints.beta.type != ConstraintType::kFree) {
                    jacobian(row, GetVariableIndex(constraints, kVarBeta)) =
                        derivs.CL_beta;
                }
                if (constraints.roll_rate.type != ConstraintType::kFree) {
                    jacobian(row, GetVariableIndex(constraints, kVarRollRate)) =
                        derivs.CL_p;
                }
                if (constraints.pitch_rate.type != ConstraintType::kFree) {
                    jacobian(row,
                             GetVariableIndex(constraints, kVarPitchRate)) =
                        derivs.CL_q;
                }
                if (constraints.yaw_rate.type != ConstraintType::kFree) {
                    jacobian(row, GetVariableIndex(constraints, kVarYawRate)) =
                        derivs.CL_r;
                }
                break;
            case ConstraintType::kPitchMoment:
                residual(row) = aero.Cm - constraints.alpha.target;
                jacobian(row, GetVariableIndex(constraints, kVarAlpha)) =
                    derivs.Cm_alpha;
                break;
            default:
                break;
        }
        ++row;
    }

    // Process beta constraint
    if (constraints.beta.type != ConstraintType::kFree) {
        switch (constraints.beta.type) {
            case ConstraintType::kFixed:
                residual(row) = state.beta - constraints.beta.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kCY:
                residual(row) = aero.CY - constraints.beta.target;
                if (constraints.alpha.type != ConstraintType::kFree) {
                    jacobian(row, GetVariableIndex(constraints, kVarAlpha)) =
                        derivs.CY_alpha;
                }
                jacobian(row, GetVariableIndex(constraints, kVarBeta)) =
                    derivs.CY_beta;
                break;
            default:
                break;
        }
        ++row;
    }

    // Process roll rate constraint
    if (constraints.roll_rate.type != ConstraintType::kFree) {
        switch (constraints.roll_rate.type) {
            case ConstraintType::kFixed:
                residual(row) = state.roll_rate - constraints.roll_rate.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kRollMoment:
                residual(row) = aero.Cl - constraints.roll_rate.target;
                jacobian(row, GetVariableIndex(constraints, kVarRollRate)) =
                    derivs.Cl_p;
                break;
            default:
                break;
        }
        ++row;
    }

    // Process pitch rate constraint
    if (constraints.pitch_rate.type != ConstraintType::kFree) {
        switch (constraints.pitch_rate.type) {
            case ConstraintType::kFixed:
                residual(row) =
                    state.pitch_rate - constraints.pitch_rate.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kPitchMoment:
                residual(row) = aero.Cm - constraints.pitch_rate.target;
                jacobian(row, GetVariableIndex(constraints, kVarPitchRate)) =
                    derivs.Cm_q;
                break;
            default:
                break;
        }
        ++row;
    }

    // Process yaw rate constraint
    if (constraints.yaw_rate.type != ConstraintType::kFree) {
        switch (constraints.yaw_rate.type) {
            case ConstraintType::kFixed:
                residual(row) = state.yaw_rate - constraints.yaw_rate.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kYawMoment:
                residual(row) = aero.Cn - constraints.yaw_rate.target;
                jacobian(row, GetVariableIndex(constraints, kVarYawRate)) =
                    derivs.Cn_r;
                break;
            default:
                break;
        }
        ++row;
    }

    // Process control surface constraints
    for (const auto& [name, constraint] : constraints.control_constraints) {
        if (constraint.type == ConstraintType::kFree)
            continue;

        // Find current deflection
        Scalar deflection = 0.0;
        for (const auto& [ctrl_name, ctrl_defl] : state.control_deflections) {
            if (ctrl_name == name) {
                deflection = ctrl_defl;
                break;
            }
        }

        switch (constraint.type) {
            case ConstraintType::kFixed:
                residual(row) = deflection - constraint.target;
                jacobian(row, row) = 1.0;
                break;
            case ConstraintType::kPitchMoment:
                residual(row) = aero.Cm - constraint.target;
                // Would need control derivatives here
                jacobian(row, row) = 1.0;  // Simplified
                break;
            case ConstraintType::kRollMoment:
                residual(row) = aero.Cl - constraint.target;
                jacobian(row, row) = 1.0;  // Simplified
                break;
            case ConstraintType::kYawMoment:
                residual(row) = aero.Cn - constraint.target;
                jacobian(row, row) = 1.0;  // Simplified
                break;
            default:
                break;
        }
        ++row;
    }
}

bool TrimSolver::UpdateState(OperatingVariables& state, const VecX& delta,
                             const TrimConstraints& constraints) const {
    Index idx = 0;

    // Newton update: state_new = state_old + delta
    // where delta = J^{-1} * (-residual)

    // Update alpha
    if (constraints.alpha.type != ConstraintType::kFree) {
        state.alpha += delta(idx);
        ++idx;
    }

    // Update beta
    if (constraints.beta.type != ConstraintType::kFree) {
        state.beta += delta(idx);
        ++idx;
    }

    // Update roll rate
    if (constraints.roll_rate.type != ConstraintType::kFree) {
        state.roll_rate += delta(idx);
        ++idx;
    }

    // Update pitch rate
    if (constraints.pitch_rate.type != ConstraintType::kFree) {
        state.pitch_rate += delta(idx);
        ++idx;
    }

    // Update yaw rate
    if (constraints.yaw_rate.type != ConstraintType::kFree) {
        state.yaw_rate += delta(idx);
        ++idx;
    }

    // Update control surfaces
    for (const auto& [name, constraint] : constraints.control_constraints) {
        if (constraint.type == ConstraintType::kFree)
            continue;

        // Find or add control deflection
        bool found = false;
        for (auto& [ctrl_name, ctrl_defl] : state.control_deflections) {
            if (ctrl_name == name) {
                ctrl_defl += delta(idx);
                found = true;
                break;
            }
        }
        if (!found) {
            state.control_deflections.emplace_back(name, delta(idx));
        }
        ++idx;
    }

    return true;
}

bool TrimSolver::CheckBounds(const OperatingVariables& state) const {
    if (std::abs(state.alpha) > options_.max_alpha)
        return false;
    if (std::abs(state.beta) > options_.max_beta)
        return false;
    if (std::abs(state.roll_rate) > options_.max_rate)
        return false;
    if (std::abs(state.pitch_rate) > options_.max_rate)
        return false;
    if (std::abs(state.yaw_rate) > options_.max_rate)
        return false;
    return true;
}

Index TrimSolver::GetVariableIndex(
    const TrimConstraints& constraints, int var_type,
    [[maybe_unused]] const std::string& control_name) const {
    Index idx = 0;

    if (var_type == kVarAlpha) {
        return idx;
    }
    if (constraints.alpha.type != ConstraintType::kFree)
        ++idx;

    if (var_type == kVarBeta) {
        return idx;
    }
    if (constraints.beta.type != ConstraintType::kFree)
        ++idx;

    if (var_type == kVarRollRate) {
        return idx;
    }
    if (constraints.roll_rate.type != ConstraintType::kFree)
        ++idx;

    if (var_type == kVarPitchRate) {
        return idx;
    }
    if (constraints.pitch_rate.type != ConstraintType::kFree)
        ++idx;

    if (var_type == kVarYawRate) {
        return idx;
    }
    if (constraints.yaw_rate.type != ConstraintType::kFree)
        ++idx;

    // Control surface indices
    if (var_type >= kVarControlBase) {
        // Would need to search for control by name
        return idx;
    }

    return idx;
}

}  // namespace pegasus
