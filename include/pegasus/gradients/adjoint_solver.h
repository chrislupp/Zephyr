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

#pragma once

/// @file adjoint_solver.h
/// @brief Adjoint solver for efficient gradient computation via implicit
/// differentiation.
///
/// The VLM problem involves an implicit linear system: A(θ)·Γ = b(θ)
/// where θ are design variables, A is the AIC matrix, and Γ is circulation.
///
/// For an objective J(Γ, θ), the total derivative is:
///   dJ/dθ = ∂J/∂θ + (∂J/∂Γ)·(dΓ/dθ)
///
/// Using the adjoint method:
///   1. Solve adjoint system: Aᵀ·λ = (∂J/∂Γ)ᵀ
///   2. Total gradient: dJ/dθ = ∂J/∂θ - λᵀ·(∂A/∂θ)·Γ + λᵀ·(∂b/∂θ)
///
/// This requires only ONE additional linear solve regardless of the number
/// of design variables, making it efficient for optimization.

#include "pegasus/core/types.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/solver/vortex_lattice_solver.h"

namespace pegasus::gradients {

/// @brief Adjoint solver for implicit differentiation of VLM linear system.
///
/// This class provides the machinery to compute gradients efficiently using
/// the adjoint method. It reuses the LU factorization from the forward solve.
///
/// Example usage:
/// @code
/// // After forward solve
/// VortexLatticeSolver solver;
/// solver.SetAircraft(aircraft);
/// auto solution = solver.Solve(condition);
///
/// // Compute force gradients (∂J/∂Γ)
/// VecX dJ_dGamma = compute_force_gradients(solution, ...);
///
/// // Adjoint solve
/// AdjointSolver adjoint(solver);
/// VecX lambda = adjoint.SolveAdjoint(dJ_dGamma);
///
/// // Now lambda can be used to compute total gradients
/// @endcode
class AdjointSolver {
public:
    /// @brief Construct adjoint solver from a primal solver.
    ///
    /// The primal solver must have a valid AIC matrix (i.e., SetAircraft
    /// and Solve must have been called).
    ///
    /// @param solver Reference to the primal VLM solver
    explicit AdjointSolver(const VortexLatticeSolver& solver);

    /// @brief Solve the adjoint system Aᵀ·λ = rhs.
    ///
    /// Uses the LU factorization from the forward solve for efficiency.
    ///
    /// @param rhs Right-hand side vector (typically ∂J/∂Γ)
    /// @return Adjoint variable λ
    [[nodiscard]] VecX SolveAdjoint(const VecX& rhs) const;

    /// @brief Solve multiple adjoint systems at once.
    ///
    /// Solves Aᵀ·Λ = RHS where each column is a separate adjoint problem.
    /// Useful when computing gradients for multiple objectives.
    ///
    /// @param rhs Matrix where each column is a right-hand side
    /// @return Matrix where each column is the corresponding adjoint solution
    [[nodiscard]] MatX SolveAdjointMultiple(const MatX& rhs) const;

    /// @brief Check if the adjoint solver is valid (has valid AIC).
    [[nodiscard]] bool IsValid() const { return is_valid_; }

    /// @brief Get the number of circulation unknowns.
    [[nodiscard]] Index GetNumUnknowns() const { return num_unknowns_; }

private:
    const VortexLatticeSolver* solver_;
    bool is_valid_{false};
    Index num_unknowns_{0};
};

/// @brief Compute gradients of aerodynamic forces with respect to circulation.
///
/// This computes ∂J/∂Γ for various objectives J (CL, CD, Cm, etc.).
/// These are the "seeds" for the adjoint solve.
///
/// @param mesh The vortex lattice mesh
/// @param solution The forward solution (circulation values)
/// @param condition Flight condition
/// @param reference_area Reference area for normalization
/// @param reference_chord Reference chord for moment normalization
/// @return Structure containing gradients for each output
struct ForceGradients {
    VecX dCL_dGamma;  ///< Gradient of lift coefficient
    VecX dCD_dGamma;  ///< Gradient of drag coefficient
    VecX dCY_dGamma;  ///< Gradient of side force coefficient
    VecX dCl_dGamma;  ///< Gradient of rolling moment coefficient
    VecX dCm_dGamma;  ///< Gradient of pitching moment coefficient
    VecX dCn_dGamma;  ///< Gradient of yawing moment coefficient
};

/// @brief Compute force gradients with respect to circulation.
///
/// Uses analytical derivatives of the Kutta-Joukowski force kernel.
///
/// @param mesh Vortex lattice mesh
/// @param solution Forward solution state
/// @param condition Flight condition
/// @param aircraft Aircraft geometry (for reference quantities)
/// @return Gradients of all force coefficients w.r.t. circulation
ForceGradients ComputeForceGradients(const VortexLattice& mesh,
                                     const SolutionState& solution,
                                     const FlightCondition& condition,
                                     const Aircraft& aircraft);

/// @brief Compute force gradients using reference quantities directly.
ForceGradients ComputeForceGradients(const VortexLattice& mesh,
                                     const SolutionState& solution,
                                     const FlightCondition& condition,
                                     const ReferenceQuantities& ref);

/// @brief Compute total derivative of objective w.r.t. flight condition.
///
/// Uses the adjoint method:
///   dJ/dα = ∂J/∂α|_Γ + λᵀ·(∂b/∂α)
///
/// where λ is the adjoint solution and ∂b/∂α is the RHS sensitivity.
///
/// @param objective Which coefficient to differentiate (0=CL, 1=CD, etc.)
/// @param solver Primal VLM solver (with solution)
/// @param condition Flight condition
/// @param aircraft Aircraft geometry
/// @return Total derivative dJ/dα (or dJ/dβ, etc.)
double ComputeTotalDerivative(int objective, const VortexLatticeSolver& solver,
                              const FlightCondition& condition,
                              const Aircraft& aircraft);

}  // namespace pegasus::gradients
