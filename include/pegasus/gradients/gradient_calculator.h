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

#pragma once

/// @file gradient_calculator.h
/// @brief Gradient computation for aerodynamic optimization using analytical
/// derivatives.
///
/// This module provides analytical differentiation of VLM aerodynamic outputs
/// with respect to flight conditions and geometry parameters, enabling
/// gradient-based optimization workflows via the adjoint method.

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/gradients/mesh_node_sensitivity.h"
#include "pegasus/gradients/output_selection.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <functional>

namespace pegasus::gradients {

/// @brief Flight condition gradient structure.
///
/// Contains gradients of an objective function with respect to flight
/// condition parameters.
struct FlightConditionGradients {
    Scalar d_alpha{0.0};       ///< Gradient w.r.t. angle of attack
    Scalar d_beta{0.0};        ///< Gradient w.r.t. sideslip angle
    Scalar d_roll_rate{0.0};   ///< Gradient w.r.t. roll rate
    Scalar d_pitch_rate{0.0};  ///< Gradient w.r.t. pitch rate
    Scalar d_yaw_rate{0.0};    ///< Gradient w.r.t. yaw rate
};

/// @brief Aerodynamic coefficient gradients with respect to a single input.
///
/// Contains sensitivities of all aerodynamic coefficients with respect
/// to one input parameter (e.g., dCL/dalpha, dCD/dalpha, etc.).
struct AeroGradients {
    Scalar dCL{0.0};   ///< Lift coefficient sensitivity
    Scalar dCD{0.0};   ///< Drag coefficient sensitivity
    Scalar dCY{0.0};   ///< Side force coefficient sensitivity
    Scalar dCl{0.0};   ///< Rolling moment coefficient sensitivity
    Scalar dCm{0.0};   ///< Pitching moment coefficient sensitivity
    Scalar dCn{0.0};   ///< Yawing moment coefficient sensitivity
    Scalar dCDi{0.0};  ///< Induced drag coefficient sensitivity
};

/// @brief Complete gradient information for all inputs.
struct FullGradients {
    AeroGradients alpha_derivs;
    AeroGradients beta_derivs;
    AeroGradients p_derivs;
    AeroGradients q_derivs;
    AeroGradients r_derivs;
    AeroGradients mach_derivs;
};

// ============================================================================
// Geometry Gradient Structures
// ============================================================================

/// @brief Geometry design variable types for sensitivity analysis.
enum class GeometryVariable {
    kChord,  ///< Chord length at section
    kTwist,  ///< Twist angle at section (degrees)
    kXle,    ///< X-coordinate of leading edge
    kYle,    ///< Y-coordinate of leading edge (spanwise position)
    kZle     ///< Z-coordinate of leading edge (dihedral)
};

/// @brief Gradients of all aero outputs w.r.t. a geometry variable at one
/// section.
///
/// These are computed using the adjoint method:
///   dJ/dθ = ∂J/∂θ|_Γ - λᵀ·(∂A/∂θ)·Γ + λᵀ·(∂b/∂θ)
///
/// where θ is a geometry parameter (chord, twist, position).
struct SectionGeometryGradients {
    std::string surface_name;
    Index section_index;

    AeroGradients chord_derivs;
    AeroGradients twist_derivs;
    AeroGradients xle_derivs;
    AeroGradients yle_derivs;
    AeroGradients zle_derivs;
};

/// @brief Complete geometry gradients for an entire aircraft.
struct GeometryGradients {
    std::vector<SectionGeometryGradients> section_gradients;

    [[nodiscard]] Index TotalDesignVariables() const {
        return section_gradients.size() * 5;
    }

    VecX GetChordGradientsFor(
        const std::function<Scalar(const AeroGradients&)>& selector) const {
        VecX result(section_gradients.size());
        for (size_t i = 0; i < section_gradients.size(); ++i) {
            result(i) = selector(section_gradients[i].chord_derivs);
        }
        return result;
    }

    VecX GetTwistGradientsFor(
        const std::function<Scalar(const AeroGradients&)>& selector) const {
        VecX result(section_gradients.size());
        for (size_t i = 0; i < section_gradients.size(); ++i) {
            result(i) = selector(section_gradients[i].twist_derivs);
        }
        return result;
    }
};

/// @brief Calculator for aerodynamic gradients using analytical derivatives.
///
/// This class provides methods to compute gradients of aerodynamic outputs
/// with respect to flight conditions and geometry parameters using the
/// adjoint method with analytical partial derivatives.
///
/// Example usage:
/// @code
/// GradientCalculator calc;
/// calc.SetAircraft(aircraft);
///
/// FlightCondition cond;
/// cond.alpha = 5.0 * kDegreesToRadians;
///
/// FullGradients grads = calc.ComputeFlightConditionDerivatives(cond);
/// // grads.alpha_derivs.dCL is dCL/dalpha
/// @endcode
class GradientCalculator {
public:
    GradientCalculator();
    ~GradientCalculator();

    /// @brief Set the aircraft configuration.
    void SetAircraft(const Aircraft& aircraft);

    /// @brief Set the reference point for moment calculations.
    void SetReferencePoint(const Vec3& ref);

    /// @brief Select which outputs to compute gradients for.
    ///
    /// Disabling unneeded outputs reduces the number of adjoint
    /// solves and related O(n^2) operations. By default all
    /// outputs are enabled.
    void SetOutputSelection(const OutputSelection& selection);

    /// @brief Compute derivatives with respect to all flight conditions.
    FullGradients ComputeFlightConditionDerivatives(
        const FlightCondition& condition);

    /// @brief Compute gradient with respect to angle of attack.
    AeroGradients ComputeAlphaDerivatives(const FlightCondition& condition);

    /// @brief Compute gradient with respect to sideslip beta.
    AeroGradients ComputeBetaDerivatives(const FlightCondition& condition);

    /// @brief Compute gradient with respect to roll rate p.
    AeroGradients ComputeRollRateDerivatives(const FlightCondition& condition);

    /// @brief Compute gradient with respect to pitch rate q.
    AeroGradients ComputePitchRateDerivatives(const FlightCondition& condition);

    /// @brief Compute gradient with respect to yaw rate r.
    AeroGradients ComputeYawRateDerivatives(const FlightCondition& condition);

    /// @brief Compute gradient with respect to Mach number.
    ///
    /// Mach affects the Prandtl-Glauert compressibility correction in the
    /// influence matrix. This uses a two-solve approach (re-solve at perturbed
    /// Mach numbers) since Mach changes the AIC matrix structure.
    AeroGradients ComputeMachDerivatives(const FlightCondition& condition);

    // ==========================================================================
    // Geometry Gradients
    // ==========================================================================

    /// @brief Compute gradients with respect to all geometry parameters.
    ///
    /// For each section in the aircraft, computes gradients of all aerodynamic
    /// coefficients with respect to chord, twist, and LE position using the
    /// adjoint method with analytical partial derivatives.
    GeometryGradients ComputeGeometryGradients(
        const FlightCondition& condition);

    /// @brief Compute gradient w.r.t. a single geometry parameter at one
    /// section.
    AeroGradients ComputeSectionGradient(const FlightCondition& condition,
                                         Index surface_index,
                                         Index section_index,
                                         GeometryVariable variable);

    // ==========================================================================
    // Direct Mesh Input
    // ==========================================================================

    /// @brief Set direct mesh input (alternative to SetAircraft).
    void SetMesh(const MeshInput& input);

    /// @brief Compute gradients of all aero outputs w.r.t. mesh node
    /// coordinates.
    ///
    /// Returns one SurfaceMeshNodeGradients per surface, containing gradients
    /// of CL, CD, CY, Cl, Cm, Cn w.r.t. each node's x,y,z coordinates.
    std::vector<SurfaceMeshNodeGradients> ComputeMeshNodeGradients(
        const FlightCondition& condition);

    /// @brief No-op (analytical derivatives don't use finite differences).
    void SetFiniteDiffStep(Scalar step);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pegasus::gradients
