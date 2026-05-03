/*
    Project Pegasus
    Version 0.1.0

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

#include "pegasus/gradients/adjoint_solver.h"

#include "pegasus/geometry/mesh_input.h"
#include "pegasus/gradients/analytical_biot_savart.h"
#include "pegasus/gradients/analytical_forces.h"

#include <cmath>
#include <vector>

namespace pegasus::gradients {

// ============================================================================
// AdjointSolver Implementation
// ============================================================================

AdjointSolver::AdjointSolver(const VortexLatticeSolver& solver)
    : solver_(&solver) {
    if (solver.IsAICValid()) {
        is_valid_ = true;
        num_unknowns_ = solver.GetAICMatrix().rows();
    }
}

VecX AdjointSolver::SolveAdjoint(const VecX& rhs) const {
    if (!is_valid_) {
        return VecX::Zero(rhs.size());
    }

    const MatX& A = solver_->GetAICMatrix();
    Eigen::PartialPivLU<MatX> lu_transpose(A.transpose());
    return lu_transpose.solve(rhs);
}

MatX AdjointSolver::SolveAdjointMultiple(const MatX& rhs) const {
    if (!is_valid_) {
        return MatX::Zero(rhs.rows(), rhs.cols());
    }

    const MatX& A = solver_->GetAICMatrix();
    Eigen::PartialPivLU<MatX> lu_transpose(A.transpose());
    return lu_transpose.solve(rhs);
}

// ============================================================================
// Force Gradient Computation
// ============================================================================

/// Extract mesh data into raw arrays for analytical kernels
void ExtractMeshData(const VortexLattice& mesh, const SolutionState& solution,
                     const Vec3& moment_ref, std::vector<double>& bound_left,
                     std::vector<double>& bound_right,
                     std::vector<double>& control_points,
                     std::vector<double>& moment_arm,
                     std::vector<double>& induced_velocity) {
    const Index n = mesh.GetNumVortices();
    bound_left.resize(3 * n);
    bound_right.resize(3 * n);
    control_points.resize(3 * n);
    moment_arm.resize(3 * n);
    induced_velocity.resize(3 * n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        bound_left[3 * i + 0] = vor.bound_vortex_left.x();
        bound_left[3 * i + 1] = vor.bound_vortex_left.y();
        bound_left[3 * i + 2] = vor.bound_vortex_left.z();

        bound_right[3 * i + 0] = vor.bound_vortex_right.x();
        bound_right[3 * i + 1] = vor.bound_vortex_right.y();
        bound_right[3 * i + 2] = vor.bound_vortex_right.z();

        control_points[3 * i + 0] = vor.control_point.x();
        control_points[3 * i + 1] = vor.control_point.y();
        control_points[3 * i + 2] = vor.control_point.z();

        moment_arm[3 * i + 0] = vor.control_point.x() - moment_ref.x();
        moment_arm[3 * i + 1] = vor.control_point.y() - moment_ref.y();
        moment_arm[3 * i + 2] = vor.control_point.z() - moment_ref.z();

        if (solution.induced_velocity.rows() > i) {
            induced_velocity[3 * i + 0] = solution.induced_velocity(i, 0);
            induced_velocity[3 * i + 1] = solution.induced_velocity(i, 1);
            induced_velocity[3 * i + 2] = solution.induced_velocity(i, 2);
        } else {
            induced_velocity[3 * i + 0] = 0.0;
            induced_velocity[3 * i + 1] = 0.0;
            induced_velocity[3 * i + 2] = 0.0;
        }
    }
}

/// @brief Compute the velocity kernel transpose applied to dJ/dV.
///
/// This computes: result[j] = sum_i K[i,j] · dJ_dV[i]
/// where K[i,j] is the velocity induced at bound vortex center i by unit
/// circulation at vortex j.
///
/// This term is needed for the correct adjoint RHS:
/// dJ/dGamma_total = dJ/dGamma|_V + K^T * dJ/dV
void ApplyVelocityKernelTranspose(const VortexLattice& mesh,
                                  const FlightCondition& condition,
                                  const std::vector<double>& dJ_dV,
                                  VecX& result) {
    const Index n = mesh.GetNumVortices();
    result.resize(n);
    result.setZero();

    Scalar mach = condition.mach;
    Scalar prandtl_glauert = (mach < 1.0) ? std::sqrt(1.0 - mach * mach) : 1.0;

    for (Index j = 0; j < n; ++j) {
        const HorseshoeVortex& vor_j = mesh.GetVortex(j);

        double sum = 0.0;

        Scalar panel_width_j =
            (vor_j.bound_vortex_right - vor_j.bound_vortex_left).norm();

        for (Index i = 0; i < n; ++i) {
            const HorseshoeVortex& vor_i = mesh.GetVortex(i);

            Vec3 field_point = vor_i.control_point;

            // Match solver's core radius logic:
            // same_component → tiny core (0.0001 * width)
            // different_component → larger core
            bool same_component =
                (vor_i.component_index == vor_j.component_index);
            Scalar core;
            if (same_component) {
                core = 0.0001 * panel_width_j;
            } else {
                Scalar chord_core =
                    mesh.GetVortexCoreRadiusChordFraction() * vor_j.local_chord;
                Scalar width_core =
                    mesh.GetVortexCoreRadiusWidthFraction() * panel_width_j;
                core = std::max(chord_core, width_core);
            }

            bool include_bound = (i != j);

            double fp[3] = {field_point.x(), field_point.y(), field_point.z()};
            double bl[3] = {vor_j.bound_vortex_left.x(),
                            vor_j.bound_vortex_left.y(),
                            vor_j.bound_vortex_left.z()};
            double br[3] = {vor_j.bound_vortex_right.x(),
                            vor_j.bound_vortex_right.y(),
                            vor_j.bound_vortex_right.z()};

            double vel[3];
            compute_horseshoe_velocity_kernel(fp, bl, br, prandtl_glauert, core,
                                              include_bound, vel);

            sum += vel[0] * dJ_dV[3 * i + 0] + vel[1] * dJ_dV[3 * i + 1] +
                   vel[2] * dJ_dV[3 * i + 2];
        }

        result(j) = sum;
    }
}

/// Core implementation of force gradient computation
static ForceGradients ComputeForceGradientsImpl(
    const VortexLattice& mesh, const SolutionState& solution,
    const FlightCondition& condition, const Vec3& moment_ref,
    double reference_area) {
    const Index n = solution.circulation.size();

    ForceGradients grads;
    grads.dCL_dGamma = VecX::Zero(n);
    grads.dCD_dGamma = VecX::Zero(n);
    grads.dCY_dGamma = VecX::Zero(n);
    grads.dCl_dGamma = VecX::Zero(n);
    grads.dCm_dGamma = VecX::Zero(n);
    grads.dCn_dGamma = VecX::Zero(n);

    if (!solution.is_valid || n == 0) {
        return grads;
    }

    // Extract mesh data into raw arrays
    std::vector<double> bound_left, bound_right, control_points, moment_arm,
        induced_velocity;
    ExtractMeshData(mesh, solution, moment_ref, bound_left, bound_right,
                    control_points, moment_arm, induced_velocity);

    // Copy circulation to contiguous array
    std::vector<double> circulation(n);
    for (Index i = 0; i < n; ++i) {
        circulation[i] = solution.circulation(i);
    }

    // The TOTAL circulation gradient is:
    //   dJ/dGamma_total = dJ/dGamma|_V + K^T * dJ/dV
    // where K is the velocity kernel at bound vortex centers.

    std::vector<double> d_circulation(n);
    std::vector<double> d_induced_velocity(3 * n);
    VecX velocity_kernel_term;

    // CL gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 0, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 0, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCL_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    // CD gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 1, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 1, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCD_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    // CY gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 2, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 2, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCY_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    // Cl gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 3, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 3, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCl_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    // Cm gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 4, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 4, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCm_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    // Cn gradient
    analytical_force_circulation_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 5, d_circulation.data());
    analytical_force_induced_velocity_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        control_points.data(), moment_arm.data(), static_cast<int>(n),
        reference_area, 5, d_induced_velocity.data());
    ApplyVelocityKernelTranspose(mesh, condition, d_induced_velocity,
                                 velocity_kernel_term);
    for (Index i = 0; i < n; ++i) {
        grads.dCn_dGamma(i) = d_circulation[i] + velocity_kernel_term(i);
    }

    return grads;
}

/// Compute force gradients using Aircraft (existing API)
ForceGradients ComputeForceGradients(const VortexLattice& mesh,
                                     const SolutionState& solution,
                                     const FlightCondition& condition,
                                     const Aircraft& aircraft) {
    return ComputeForceGradientsImpl(mesh, solution, condition,
                                     aircraft.GetMomentReferencePoint(),
                                     aircraft.GetReferenceArea());
}

/// Compute force gradients using ReferenceQuantities (mesh-based API)
ForceGradients ComputeForceGradients(const VortexLattice& mesh,
                                     const SolutionState& solution,
                                     const FlightCondition& condition,
                                     const ReferenceQuantities& ref) {
    return ComputeForceGradientsImpl(mesh, solution, condition,
                                     ref.moment_reference_point,
                                     ref.reference_area);
}

// ============================================================================
// Total Derivative Computation
// ============================================================================

double ComputeTotalDerivative(int objective, const VortexLatticeSolver& solver,
                              const FlightCondition& condition,
                              const Aircraft& aircraft) {
    const SolutionState& solution = solver.GetSolution();
    const VortexLattice& mesh = solver.GetMesh();

    if (!solution.is_valid) {
        return 0.0;
    }

    // Step 1: Compute force gradients ∂J/∂Γ
    ForceGradients force_grads =
        ComputeForceGradients(mesh, solution, condition, aircraft);

    // Select the appropriate gradient based on objective
    VecX dJ_dGamma;
    switch (objective) {
        case 0:
            dJ_dGamma = force_grads.dCL_dGamma;
            break;
        case 1:
            dJ_dGamma = force_grads.dCD_dGamma;
            break;
        case 2:
            dJ_dGamma = force_grads.dCY_dGamma;
            break;
        case 3:
            dJ_dGamma = force_grads.dCl_dGamma;
            break;
        case 4:
            dJ_dGamma = force_grads.dCm_dGamma;
            break;
        case 5:
            dJ_dGamma = force_grads.dCn_dGamma;
            break;
        default:
            return 0.0;
    }

    // Step 2: Solve adjoint system Aᵀ·λ = ∂J/∂Γ
    AdjointSolver adjoint(solver);
    VecX lambda = adjoint.SolveAdjoint(dJ_dGamma);

    // Step 3: Compute ∂b/∂α analytically
    // V_inf = (cos(α)cos(β), sin(β), sin(α)cos(β))
    // ∂V_inf/∂α = (-sin(α)cos(β), 0, cos(α)cos(β))
    double dv_dalpha[3];
    compute_freestream_alpha_gradient(condition.alpha, condition.beta,
                                      dv_dalpha);
    Vec3 dV_dalpha(dv_dalpha[0], dv_dalpha[1], dv_dalpha[2]);

    const Index n = solution.circulation.size();
    VecX db_dalpha(n);
    for (Index i = 0; i < n; ++i) {
        Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
        db_dalpha(i) = -normal.dot(dV_dalpha);
    }

    double implicit_term = lambda.dot(db_dalpha);

    return implicit_term;
}

}  // namespace pegasus::gradients
