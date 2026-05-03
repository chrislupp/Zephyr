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

#ifndef PEGASUS_AERO_STABILITY_DERIVATIVES_H_
#define PEGASUS_AERO_STABILITY_DERIVATIVES_H_

#include "pegasus/core/types.h"

#include <string>
#include <vector>

namespace pegasus {

/// Control surface derivative data
struct ControlDerivatives {
    std::string name;        ///< Control surface name
    Scalar deflection{0.0};  ///< Current deflection (rad)
    Scalar gain{1.0};        ///< Control gain

    // Force derivatives w.r.t. control deflection (stability axis)
    Scalar CL_d{0.0};  ///< dCL/d(deflection)
    Scalar CY_d{0.0};  ///< dCY/d(deflection)
    Scalar CD_d{0.0};  ///< dCD/d(deflection)

    // Moment derivatives w.r.t. control deflection (stability axis)
    Scalar Cl_d{0.0};  ///< dCl/d(deflection) - roll moment
    Scalar Cm_d{0.0};  ///< dCm/d(deflection) - pitch moment
    Scalar Cn_d{0.0};  ///< dCn/d(deflection) - yaw moment

    // Force derivatives (body axis)
    Scalar CX_d{0.0};       ///< dCX/d(deflection)
    Scalar CY_d_body{0.0};  ///< dCY/d(deflection) body axis
    Scalar CZ_d{0.0};       ///< dCZ/d(deflection)

    // Hinge moment derivative
    Scalar Ch_d{0.0};  ///< dCh/d(deflection)
};

/// Stability derivatives in stability axis system
///
/// Sign conventions follow standard aerodynamic notation:
/// - α (alpha): angle of attack, positive nose up
/// - β (beta): sideslip angle, positive nose right
/// - p: roll rate, positive right wing down
/// - q: pitch rate, positive nose up
/// - r: yaw rate, positive nose right
///
/// Rate derivatives are non-dimensionalized:
/// - Roll/yaw derivatives: by (2V/b) where b is span
/// - Pitch derivatives: by (2V/c) where c is chord
struct StabilityDerivatives {
    bool is_valid{false};

    // Reference quantities used for non-dimensionalization
    Scalar reference_area{1.0};
    Scalar reference_chord{1.0};
    Scalar reference_span{1.0};
    Scalar velocity{1.0};

    // Flight condition at which derivatives were computed
    Scalar alpha{0.0};  ///< Angle of attack (rad)
    Scalar beta{0.0};   ///< Sideslip angle (rad)

    // =========================================================================
    // Angle of attack derivatives (per radian)
    // =========================================================================

    Scalar CL_alpha{0.0};  ///< Lift curve slope
    Scalar CY_alpha{0.0};  ///< Side force due to alpha
    Scalar CD_alpha{0.0};  ///< Drag due to alpha
    Scalar Cl_alpha{0.0};  ///< Roll moment due to alpha
    Scalar Cm_alpha{0.0};  ///< Pitch moment due to alpha (static stability)
    Scalar Cn_alpha{0.0};  ///< Yaw moment due to alpha

    // =========================================================================
    // Sideslip derivatives (per radian)
    // =========================================================================

    Scalar CL_beta{0.0};  ///< Lift due to sideslip
    Scalar CY_beta{0.0};  ///< Side force due to sideslip (weathercock)
    Scalar CD_beta{0.0};  ///< Drag due to sideslip
    Scalar Cl_beta{0.0};  ///< Roll due to sideslip (dihedral effect)
    Scalar Cm_beta{0.0};  ///< Pitch due to sideslip
    Scalar Cn_beta{0.0};  ///< Yaw due to sideslip (directional stability)

    // =========================================================================
    // Roll rate derivatives (non-dimensional: p * b / 2V)
    // =========================================================================

    Scalar CL_p{0.0};  ///< Lift due to roll rate
    Scalar CY_p{0.0};  ///< Side force due to roll rate
    Scalar CD_p{0.0};  ///< Drag due to roll rate
    Scalar Cl_p{0.0};  ///< Roll damping
    Scalar Cm_p{0.0};  ///< Pitch due to roll rate
    Scalar Cn_p{0.0};  ///< Yaw due to roll rate

    // =========================================================================
    // Pitch rate derivatives (non-dimensional: q * c / 2V)
    // =========================================================================

    Scalar CL_q{0.0};  ///< Lift due to pitch rate
    Scalar CY_q{0.0};  ///< Side force due to pitch rate
    Scalar CD_q{0.0};  ///< Drag due to pitch rate
    Scalar Cl_q{0.0};  ///< Roll due to pitch rate
    Scalar Cm_q{0.0};  ///< Pitch damping
    Scalar Cn_q{0.0};  ///< Yaw due to pitch rate

    // =========================================================================
    // Yaw rate derivatives (non-dimensional: r * b / 2V)
    // =========================================================================

    Scalar CL_r{0.0};  ///< Lift due to yaw rate
    Scalar CY_r{0.0};  ///< Side force due to yaw rate
    Scalar CD_r{0.0};  ///< Drag due to yaw rate
    Scalar Cl_r{0.0};  ///< Roll due to yaw rate
    Scalar Cm_r{0.0};  ///< Pitch due to yaw rate
    Scalar Cn_r{0.0};  ///< Yaw damping

    // =========================================================================
    // Control surface derivatives
    // =========================================================================

    std::vector<ControlDerivatives> control_derivatives;

    // =========================================================================
    // Neutral point and stability margins
    // =========================================================================

    Scalar neutral_point_x{0.0};  ///< X location of neutral point
    Scalar static_margin{0.0};    ///< (x_np - x_cg) / c_ref
};

/// Stability derivatives in body axis system
///
/// Body axis system:
/// - X: forward (positive out the nose)
/// - Y: right wing (positive out right wing)
/// - Z: down (positive down)
struct BodyAxisDerivatives {
    bool is_valid{false};

    // Reference quantities
    Scalar reference_area{1.0};
    Scalar reference_chord{1.0};
    Scalar reference_span{1.0};
    Scalar velocity{1.0};

    // Flight condition
    Scalar alpha{0.0};
    Scalar beta{0.0};

    // =========================================================================
    // Angle of attack derivatives (per radian)
    // =========================================================================

    Scalar CX_alpha{0.0};  ///< Axial force due to alpha
    Scalar CY_alpha{0.0};  ///< Side force due to alpha
    Scalar CZ_alpha{0.0};  ///< Normal force due to alpha
    Scalar Cl_alpha{0.0};  ///< Roll moment due to alpha
    Scalar Cm_alpha{0.0};  ///< Pitch moment due to alpha
    Scalar Cn_alpha{0.0};  ///< Yaw moment due to alpha

    // =========================================================================
    // Sideslip derivatives (per radian)
    // =========================================================================

    Scalar CX_beta{0.0};  ///< Axial force due to beta
    Scalar CY_beta{0.0};  ///< Side force due to beta
    Scalar CZ_beta{0.0};  ///< Normal force due to beta
    Scalar Cl_beta{0.0};  ///< Roll moment due to beta
    Scalar Cm_beta{0.0};  ///< Pitch moment due to beta
    Scalar Cn_beta{0.0};  ///< Yaw moment due to beta

    // =========================================================================
    // Roll rate derivatives (non-dimensional: p * b / 2V)
    // =========================================================================

    Scalar CX_p{0.0};
    Scalar CY_p{0.0};
    Scalar CZ_p{0.0};
    Scalar Cl_p{0.0};
    Scalar Cm_p{0.0};
    Scalar Cn_p{0.0};

    // =========================================================================
    // Pitch rate derivatives (non-dimensional: q * c / 2V)
    // =========================================================================

    Scalar CX_q{0.0};
    Scalar CY_q{0.0};
    Scalar CZ_q{0.0};
    Scalar Cl_q{0.0};
    Scalar Cm_q{0.0};
    Scalar Cn_q{0.0};

    // =========================================================================
    // Yaw rate derivatives (non-dimensional: r * b / 2V)
    // =========================================================================

    Scalar CX_r{0.0};
    Scalar CY_r{0.0};
    Scalar CZ_r{0.0};
    Scalar Cl_r{0.0};
    Scalar Cm_r{0.0};
    Scalar Cn_r{0.0};
};

// Forward declarations
class Aircraft;
class VortexLatticeSolver;
struct FlightCondition;

/// Calculator class for computing stability derivatives from VLM solution.
///
/// Uses the unit circulation responses from the solver to efficiently
/// compute derivatives without re-solving the system. The derivatives
/// are computed by applying the chain rule to transform from velocity
/// component derivatives to angle derivatives.
class StabilityDerivativeCalculator {
public:
    StabilityDerivativeCalculator() = default;

    /// Compute stability axis derivatives at a given flight condition.
    ///
    /// @param aircraft Aircraft geometry (for reference quantities)
    /// @param solver VLM solver with unit responses computed
    /// @param condition Flight condition (defines operating point)
    /// @return Stability derivatives in stability axis system
    StabilityDerivatives ComputeDerivatives(const Aircraft& aircraft,
                                            VortexLatticeSolver& solver,
                                            const FlightCondition& condition);

    /// Compute body axis derivatives at a given flight condition.
    ///
    /// @param aircraft Aircraft geometry (for reference quantities)
    /// @param solver VLM solver with unit responses computed
    /// @param condition Flight condition (defines operating point)
    /// @return Stability derivatives in body axis system
    BodyAxisDerivatives ComputeBodyAxisDerivatives(
        const Aircraft& aircraft, VortexLatticeSolver& solver,
        const FlightCondition& condition);
};

}  // namespace pegasus

#endif  // PEGASUS_AERO_STABILITY_DERIVATIVES_H_
