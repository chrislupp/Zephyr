/*
    Project VLM
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

#ifndef VLM_SOLVER_INFLUENCE_MATRIX_H_
#define VLM_SOLVER_INFLUENCE_MATRIX_H_

#include "vlm/core/types.h"
#include "vlm/solver/vortex_lattice.h"

namespace vlm {

/// @brief Computes induced velocity from a horseshoe vortex using Biot-Savart
/// law.
///
/// Uses the Leishman R^4 core model for regularization near the vortex line.
/// Includes Prandtl-Glauert compressibility correction.
///
/// @param field_point Point where velocity is computed
/// @param bound_left Left endpoint of bound vortex
/// @param bound_right Right endpoint of bound vortex
/// @param prandtl_glauert_factor sqrt(1 - M^2), compression factor
/// @param core_radius Finite core radius for regularization
/// @param include_bound Whether to include the bound vortex contribution
/// @return Induced velocity vector (per unit circulation)
Vec3 ComputeHorseshoeVortexVelocity(const Vec3& field_point,
                                    const Vec3& bound_left,
                                    const Vec3& bound_right,
                                    Scalar prandtl_glauert_factor,
                                    Scalar core_radius,
                                    bool include_bound = true);

/// @brief Configuration for symmetry planes in AIC computation.
struct SymmetryConfig {
    /// Y-plane symmetry type (0=none, 1=symmetric, -1=antisymmetric)
    SymmetryType y_symmetry{SymmetryType::kNone};

    /// Y coordinate of symmetry plane
    Scalar y_symmetry_plane{0.0};

    /// Z-plane symmetry type (0=none, 1=symmetric, -1=antisymmetric)
    SymmetryType z_symmetry{SymmetryType::kNone};

    /// Z coordinate of symmetry plane
    Scalar z_symmetry_plane{0.0};
};

/// @brief Builds aerodynamic influence coefficient matrices.
///
/// Computes the velocity influence matrix (WC_GAM) and normalwash influence
/// matrix (AIC) for a vortex lattice mesh. The AIC matrix relates vortex
/// strengths to normal wash at control points.
///
/// Reference: Katz & Plotkin, "Low-Speed Aerodynamics", Chapter 12
class InfluenceMatrixBuilder {
public:
    /// @brief Constructs builder with mesh reference.
    explicit InfluenceMatrixBuilder(const VortexLattice& mesh);

    /// @name Configuration
    /// @{

    /// @brief Sets the Mach number for compressibility correction.
    void SetMachNumber(Scalar mach);

    /// @brief Gets the current Mach number.
    Scalar GetMachNumber() const { return mach_number_; }

    /// @brief Gets the Prandtl-Glauert compressibility factor.
    Scalar GetPrandtlGlauertFactor() const { return prandtl_glauert_factor_; }

    /// @brief Sets symmetry configuration.
    void SetSymmetry(const SymmetryConfig& config) { symmetry_ = config; }

    /// @brief Gets the symmetry configuration.
    const SymmetryConfig& GetSymmetry() const { return symmetry_; }

    /// @brief Sets vortex core radius as fraction of local chord.
    void SetCoreRadiusChordFraction(Scalar fraction) {
        core_radius_chord_fraction_ = fraction;
    }

    /// @brief Sets vortex core radius as fraction of panel width.
    void SetCoreRadiusWidthFraction(Scalar fraction) {
        core_radius_width_fraction_ = fraction;
    }

    /// @}

    /// @name Matrix Building
    /// @{

    /// @brief Builds the velocity influence matrix.
    ///
    /// Computes a 3N x N matrix where element (3i+k, j) gives the k-th
    /// velocity component at control point i due to unit circulation on
    /// horseshoe vortex j.
    ///
    /// @return Velocity influence matrix (3*num_vortices x num_vortices)
    MatX BuildVelocityInfluenceMatrix() const;

    /// @brief Builds the normalwash (AIC) influence matrix.
    ///
    /// Projects velocity influence onto panel normals to get the N x N
    /// aerodynamic influence coefficient matrix.
    ///
    /// @return AIC matrix (num_vortices x num_vortices)
    MatX BuildNormalwashMatrix() const;

    /// @brief Builds AIC matrix with optional velocity matrix output.
    ///
    /// @param[out] velocity_matrix If not null, stores velocity influence
    /// matrix
    /// @return AIC matrix
    MatX BuildNormalwashMatrix(MatX* velocity_matrix) const;

    /// @}

    /// @name Velocity Computation
    /// @{

    /// @brief Computes induced velocity at a point from all vortices.
    ///
    /// @param field_point Point where velocity is computed
    /// @param circulation Vector of circulation strengths
    /// @return Induced velocity vector
    Vec3 ComputeInducedVelocity(const Vec3& field_point,
                                const VecX& circulation) const;

    /// @brief Computes induced velocity from a single horseshoe vortex.
    ///
    /// @param field_point Point where velocity is computed
    /// @param vortex_index Index of the horseshoe vortex
    /// @param include_bound Whether to include bound vortex contribution
    /// @return Induced velocity vector (per unit circulation)
    Vec3 ComputeSingleVortexVelocity(const Vec3& field_point,
                                     Index vortex_index,
                                     bool include_bound = true) const;

    /// @}

private:
    /// @brief Computes core radius for a vortex element.
    Scalar ComputeCoreRadius(Index vortex_index, bool same_component) const;

    /// @brief Computes velocity from vortex with Y-symmetry image.
    Vec3 ComputeVelocityWithYImage(const Vec3& field_point, Index vortex_index,
                                   bool include_bound,
                                   Scalar core_radius) const;

    /// @brief Computes velocity from vortex with Z-symmetry image.
    Vec3 ComputeVelocityWithZImage(const Vec3& field_point, Index vortex_index,
                                   bool include_bound,
                                   Scalar core_radius) const;

    const VortexLattice& mesh_;
    Scalar mach_number_{0.0};
    Scalar prandtl_glauert_factor_{1.0};
    SymmetryConfig symmetry_;
    Scalar core_radius_chord_fraction_{0.0001};
    Scalar core_radius_width_fraction_{0.0001};
};

}  // namespace vlm

#endif  // VLM_SOLVER_INFLUENCE_MATRIX_H_
