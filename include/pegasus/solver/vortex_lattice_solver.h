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

#ifndef PEGASUS_SOLVER_VORTEX_LATTICE_SOLVER_H_
#define PEGASUS_SOLVER_VORTEX_LATTICE_SOLVER_H_

#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/solver/influence_matrix.h"
#include "pegasus/solver/vortex_lattice.h"

#include <memory>
#include <string>
#include <vector>

namespace pegasus {

/// @brief Flight condition for VLM analysis.
struct FlightCondition {
    /// Angle of attack in radians
    Scalar alpha{0.0};

    /// Sideslip angle in radians
    Scalar beta{0.0};

    /// Roll rate (p) normalized by 2V/b
    Scalar roll_rate{0.0};

    /// Pitch rate (q) normalized by 2V/c
    Scalar pitch_rate{0.0};

    /// Yaw rate (r) normalized by 2V/b
    Scalar yaw_rate{0.0};

    /// Freestream Mach number
    Scalar mach{0.0};

    /// Dynamic pressure (0.5 * rho * V^2)
    Scalar dynamic_pressure{1.0};

    /// Control surface deflections (radians), indexed by control name
    std::vector<std::pair<std::string, Scalar>> control_deflections;
};

/// @brief Results from VLM solution.
struct SolutionState {
    /// Vortex circulation strengths
    VecX circulation;

    /// Circulation per unit freestream velocity (u component)
    VecX circulation_per_u;

    /// Circulation per unit freestream velocity (v component)
    VecX circulation_per_v;

    /// Circulation per unit freestream velocity (w component)
    VecX circulation_per_w;

    /// Circulation per unit roll rate
    VecX circulation_per_p;

    /// Circulation per unit pitch rate
    VecX circulation_per_q;

    /// Circulation per unit yaw rate
    VecX circulation_per_r;

    /// Induced velocity at control points
    MatX induced_velocity;

    /// Whether solution is valid
    bool is_valid{false};
};

/// @brief Main vortex lattice method solver.
///
/// Solves for circulation distribution on a discretized lifting surface
/// using the vortex lattice method. The solver:
/// 1. Generates mesh from aircraft geometry
/// 2. Builds aerodynamic influence coefficient (AIC) matrix
/// 3. Factors AIC matrix (LU decomposition)
/// 4. Solves for unit freestream responses
/// 5. Combines responses for actual flight condition
///
/// The solver caches the AIC factorization for efficient analysis
/// at multiple flight conditions (same geometry, different alpha/beta).
class VortexLatticeSolver {
public:
    VortexLatticeSolver();
    ~VortexLatticeSolver();

    /// @name Setup
    /// @{

    /// @brief Sets the aircraft configuration.
    ///
    /// Generates the vortex lattice mesh from the aircraft geometry.
    /// Invalidates any cached AIC matrix.
    void SetAircraft(const Aircraft& aircraft);

    /// @brief Sets direct mesh input (alternative to SetAircraft).
    ///
    /// Generates the vortex lattice from structured mesh node grids.
    /// Reference quantities are taken from the MeshInput.
    void SetMesh(const MeshInput& input);

    /// @brief Gets reference quantities (set by SetMesh or extracted from
    /// Aircraft).
    [[nodiscard]] const ReferenceQuantities& GetReferenceQuantities() const {
        return ref_quantities_;
    }

    /// @brief Gets the generated mesh.
    [[nodiscard]] const VortexLattice& GetMesh() const { return mesh_; }

    /// @brief Gets a mutable reference to the mesh.
    VortexLattice& GetMeshMutable() { return mesh_; }

    /// @brief Sets symmetry configuration.
    void SetSymmetry(const SymmetryConfig& config);

    /// @brief Gets the symmetry configuration.
    [[nodiscard]] const SymmetryConfig& GetSymmetry() const {
        return symmetry_;
    }

    /// @brief Sets the reference point for moment calculations.
    void SetReferencePoint(const Vec3& ref) { reference_point_ = ref; }

    /// @brief Gets the reference point.
    [[nodiscard]] const Vec3& GetReferencePoint() const {
        return reference_point_;
    }

    /// @}

    /// @name Solver Control
    /// @{

    /// @brief Rebuilds the AIC matrix and factorization.
    ///
    /// Call this after modifying the mesh or changing Mach number.
    void RebuildAIC();

    /// @brief Solves for unit freestream responses.
    ///
    /// Computes circulation distribution for:
    /// - Unit freestream velocity in x, y, z
    /// - Unit rotation rates p, q, r
    ///
    /// These are used to efficiently compute solutions for any
    /// flight condition without re-solving the system.
    void SolveUnitResponses();

    /// @brief Checks if AIC matrix is valid.
    [[nodiscard]] bool IsAICValid() const { return aic_valid_; }

    /// @brief Checks if unit responses are computed.
    [[nodiscard]] bool AreUnitResponsesValid() const {
        return unit_responses_valid_;
    }

    /// @}

    /// @name Solution
    /// @{

    /// @brief Solves for circulation at given flight condition.
    ///
    /// @param condition Flight condition (alpha, beta, rates, etc.)
    /// @return Solution state with circulation distribution
    SolutionState Solve(const FlightCondition& condition);

    /// @brief Gets the current solution state.
    [[nodiscard]] const SolutionState& GetSolution() const { return solution_; }

    /// @brief Computes freestream velocity vector from flight condition.
    ///
    /// @param condition Flight condition
    /// @return Unit freestream velocity in body axes
    static Vec3 ComputeFreestreamDirection(const FlightCondition& condition);

    /// @}

    /// @name Matrix Access
    /// @{

    /// @brief Gets the AIC matrix.
    [[nodiscard]] const MatX& GetAICMatrix() const { return aic_matrix_; }

    /// @brief Gets the velocity influence matrix.
    [[nodiscard]] const MatX& GetVelocityInfluenceMatrix() const {
        return velocity_influence_matrix_;
    }

    /// @}

private:
    /// @brief Builds right-hand side for unit freestream component.
    [[nodiscard]] VecX BuildUnitFreestreamRHS(int component) const;

    /// @brief Builds right-hand side for unit rotation component.
    [[nodiscard]] VecX BuildUnitRotationRHS(int component) const;

    /// @brief Combines unit responses for flight condition.
    [[nodiscard]] VecX CombineUnitResponses(
        const FlightCondition& condition) const;

    VortexLattice mesh_;
    SymmetryConfig symmetry_;
    Vec3 reference_point_{0, 0, 0};
    ReferenceQuantities ref_quantities_;

    // AIC matrix and factorization
    MatX aic_matrix_;
    MatX velocity_influence_matrix_;
    Eigen::PartialPivLU<MatX> aic_lu_;
    bool aic_valid_{false};

    // Unit freestream responses (6 columns: u, v, w, p, q, r)
    MatX unit_circulation_;
    bool unit_responses_valid_{false};

    // Current solution
    SolutionState solution_;
    Scalar current_mach_{0.0};
};

}  // namespace pegasus

#endif  // PEGASUS_SOLVER_VORTEX_LATTICE_SOLVER_H_
