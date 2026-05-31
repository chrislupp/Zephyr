/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

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

#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>
#include <stdexcept>

namespace pegasus {

VortexLatticeSolver::VortexLatticeSolver() = default;

VortexLatticeSolver::~VortexLatticeSolver() = default;

void VortexLatticeSolver::SetAircraft(const Aircraft& aircraft) {
    mesh_.GenerateFromAircraft(aircraft);

    // Use aircraft moment reference point for rotation calculations
    reference_point_ = aircraft.GetMomentReferencePoint();

    // Invalidate cached matrices
    aic_valid_ = false;
    unit_responses_valid_ = false;
}

void VortexLatticeSolver::SetMesh(const MeshInput& input) {
    mesh_.GenerateFromMesh(input);

    ref_quantities_ = input.reference;
    reference_point_ = input.reference.moment_reference_point;

    // Invalidate cached matrices
    aic_valid_ = false;
    unit_responses_valid_ = false;
}

void VortexLatticeSolver::SetSymmetry(const SymmetryConfig& config) {
    symmetry_ = config;

    // Invalidate cached matrices
    aic_valid_ = false;
    unit_responses_valid_ = false;
}

void VortexLatticeSolver::RebuildAIC() {
    if (mesh_.GetNumVortices() == 0) {
        aic_valid_ = false;
        return;
    }

    InfluenceMatrixBuilder builder(mesh_);
    builder.SetMachNumber(current_mach_);
    builder.SetSymmetry(symmetry_);

    // Build normalwash matrix and velocity influence matrix
    aic_matrix_ = builder.BuildNormalwashMatrix(&velocity_influence_matrix_);

    // LU factorization for efficient solving
    aic_lu_.compute(aic_matrix_);

    aic_valid_ = true;
    unit_responses_valid_ = false;  // Need to recompute unit responses
}

VecX VortexLatticeSolver::BuildUnitFreestreamRHS(int component) const {
    Index n = mesh_.GetNumVortices();
    VecX rhs = VecX::Zero(n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh_.GetVortex(i);

        // Only apply to vortices affected by freestream
        if (!vor.affected_by_freestream) {
            continue;
        }

        // Unit freestream velocity in specified component
        Vec3 v_unit = Vec3::Zero();
        v_unit(component) = 1.0;

        // RHS = -V_inf . n
        rhs(i) = -v_unit.dot(vor.normal_at_control_point);
    }

    return rhs;
}

VecX VortexLatticeSolver::BuildUnitRotationRHS(int component) const {
    Index n = mesh_.GetNumVortices();
    VecX rhs = VecX::Zero(n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh_.GetVortex(i);

        // Only apply to vortices affected by freestream/rotation
        if (!vor.affected_by_freestream) {
            continue;
        }

        // Unit rotation rate in specified component
        Vec3 w_unit = Vec3::Zero();
        w_unit(component) = 1.0;

        // r = control_point - reference_point
        Vec3 r = vor.control_point - reference_point_;

        // V = w x r (velocity from rotation)
        Vec3 v_rot = w_unit.cross(r);

        // RHS = -V_rot . n
        rhs(i) = -v_rot.dot(vor.normal_at_control_point);
    }

    return rhs;
}

void VortexLatticeSolver::SolveUnitResponses() {
    if (!aic_valid_) {
        RebuildAIC();
    }

    if (mesh_.GetNumVortices() == 0) {
        unit_responses_valid_ = false;
        return;
    }

    Index n = mesh_.GetNumVortices();

    // 6 unit conditions: u, v, w, p, q, r
    unit_circulation_.resize(n, 6);

    // Unit freestream velocities (u, v, w)
    for (int k = 0; k < 3; ++k) {
        VecX rhs = BuildUnitFreestreamRHS(k);
        unit_circulation_.col(k) = aic_lu_.solve(rhs);
    }

    // Unit rotation rates (p, q, r)
    for (int k = 0; k < 3; ++k) {
        VecX rhs = BuildUnitRotationRHS(k);
        unit_circulation_.col(3 + k) = aic_lu_.solve(rhs);
    }

    unit_responses_valid_ = true;
}

Vec3 VortexLatticeSolver::ComputeFreestreamDirection(
    const FlightCondition& condition) {
    // Body-axis freestream velocity direction
    // For alpha (pitch): V_inf in x-z plane, rotated by alpha about y
    // For beta (yaw): V_inf has y-component
    Scalar ca = std::cos(condition.alpha);
    Scalar sa = std::sin(condition.alpha);
    Scalar cb = std::cos(condition.beta);
    Scalar sb = std::sin(condition.beta);

    // Unit freestream direction in body axes
    return {ca * cb, sb, sa * cb};
}

VecX VortexLatticeSolver::CombineUnitResponses(
    const FlightCondition& condition) const {
    // Freestream velocity components
    Vec3 v_inf = ComputeFreestreamDirection(condition);

    // Rotation rates (normalized)
    Vec3 w_rot(condition.roll_rate, condition.pitch_rate, condition.yaw_rate);

    // Combine: Gamma = Gamma_u * u + Gamma_v * v + Gamma_w * w
    //                + Gamma_p * p + Gamma_q * q + Gamma_r * r
    VecX circulation = unit_circulation_.col(0) * v_inf.x() +
                       unit_circulation_.col(1) * v_inf.y() +
                       unit_circulation_.col(2) * v_inf.z() +
                       unit_circulation_.col(3) * w_rot.x() +
                       unit_circulation_.col(4) * w_rot.y() +
                       unit_circulation_.col(5) * w_rot.z();

    return circulation;
}

SolutionState VortexLatticeSolver::Solve(const FlightCondition& condition) {
    SolutionState state;

    if (mesh_.GetNumVortices() == 0) {
        state.is_valid = false;
        return state;
    }

    // Check if we need to rebuild AIC for new Mach number
    if (std::abs(condition.mach - current_mach_) > kEpsilon) {
        current_mach_ = condition.mach;
        aic_valid_ = false;
        unit_responses_valid_ =
            false;  // AIC changed, so responses are invalid too
    }

    // Ensure unit responses are computed
    if (!unit_responses_valid_) {
        SolveUnitResponses();
    }

    // Combine unit responses for this flight condition
    state.circulation = CombineUnitResponses(condition);

    // Store unit response components for later use
    state.circulation_per_u = unit_circulation_.col(0);
    state.circulation_per_v = unit_circulation_.col(1);
    state.circulation_per_w = unit_circulation_.col(2);
    state.circulation_per_p = unit_circulation_.col(3);
    state.circulation_per_q = unit_circulation_.col(4);
    state.circulation_per_r = unit_circulation_.col(5);

    // Compute induced velocities at control points
    Index n = mesh_.GetNumVortices();
    state.induced_velocity.resize(n, 3);

    for (Index i = 0; i < n; ++i) {
        Vec3 vel = Vec3::Zero();
        for (Index j = 0; j < n; ++j) {
            vel.x() +=
                velocity_influence_matrix_(3 * i + 0, j) * state.circulation(j);
            vel.y() +=
                velocity_influence_matrix_(3 * i + 1, j) * state.circulation(j);
            vel.z() +=
                velocity_influence_matrix_(3 * i + 2, j) * state.circulation(j);
        }
        state.induced_velocity.row(i) = vel.transpose();
    }

    state.is_valid = true;
    solution_ = state;

    return state;
}

VecX VortexLatticeSolver::SolveForward(const VecX& rhs) const {
    return aic_lu_.solve(rhs);
}

}  // namespace pegasus
