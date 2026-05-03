// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_TRIM_FLIGHT_CONSTRAINTS_H_
#define VLM_TRIM_FLIGHT_CONSTRAINTS_H_

#include "vlm/core/types.h"

#include <string>
#include <vector>

namespace vlm {

/// Type of constraint on an operating variable
enum class ConstraintType {
    kFree,         ///< Variable is free (no constraint)
    kFixed,        ///< Variable is fixed at specified value
    kCL,           ///< Constrain to achieve target CL
    kCY,           ///< Constrain to achieve target side force CY
    kRollMoment,   ///< Constrain to achieve target roll moment Cl
    kPitchMoment,  ///< Constrain to achieve target pitch moment Cm
    kYawMoment,    ///< Constrain to achieve target yaw moment Cn
};

/// Constraint on a single operating variable
struct VariableConstraint {
    ConstraintType type{ConstraintType::kFree};
    Scalar target{0.0};  ///< Target value for constraint
};

/// Operating variables that can be adjusted during trim
struct OperatingVariables {
    Scalar alpha{0.0};       ///< Angle of attack (rad)
    Scalar beta{0.0};        ///< Sideslip angle (rad)
    Scalar roll_rate{0.0};   ///< Roll rate p*b/2V (non-dimensional)
    Scalar pitch_rate{0.0};  ///< Pitch rate q*c/2V (non-dimensional)
    Scalar yaw_rate{0.0};    ///< Yaw rate r*b/2V (non-dimensional)

    /// Control surface deflections (rad), indexed by control name
    std::vector<std::pair<std::string, Scalar>> control_deflections;

    /// Get number of total variables (5 base + controls)
    Index GetNumVariables() const {
        return 5 + static_cast<Index>(control_deflections.size());
    }
};

/// Constraints on operating variables for trim solution
///
/// The trim solver adjusts free variables to satisfy constraints.
/// Each operating variable (alpha, beta, p, q, r, controls) can be:
/// - Free: adjusted by the solver
/// - Fixed: held at specified value
/// - Constrained: adjusted to achieve a target (CL, Cm, etc.)
///
/// Example usage:
/// @code
/// TrimConstraints constraints;
/// // Fix beta at zero, constrain alpha to achieve CL=0.5
/// constraints.alpha.type = ConstraintType::kCL;
/// constraints.alpha.target = 0.5;
/// constraints.beta.type = ConstraintType::kFixed;
/// constraints.beta.target = 0.0;
/// // Let rates be free (will converge to zero for level flight)
/// constraints.roll_rate.type = ConstraintType::kFixed;
/// constraints.roll_rate.target = 0.0;
/// // ...
/// @endcode
struct TrimConstraints {
    /// Constraints on base operating variables
    VariableConstraint alpha;
    VariableConstraint beta;
    VariableConstraint roll_rate;
    VariableConstraint pitch_rate;
    VariableConstraint yaw_rate;

    /// Constraints on control surfaces, indexed by control name
    std::vector<std::pair<std::string, VariableConstraint>> control_constraints;

    /// Get number of constrained variables
    Index GetNumConstraints() const {
        Index count = 0;
        if (alpha.type != ConstraintType::kFree)
            ++count;
        if (beta.type != ConstraintType::kFree)
            ++count;
        if (roll_rate.type != ConstraintType::kFree)
            ++count;
        if (pitch_rate.type != ConstraintType::kFree)
            ++count;
        if (yaw_rate.type != ConstraintType::kFree)
            ++count;
        for (const auto& cc : control_constraints) {
            if (cc.second.type != ConstraintType::kFree)
                ++count;
        }
        return count;
    }
};

/// Predefined trim configurations for common flight conditions
namespace TrimPresets {

/// Level flight: constrain CL, fix beta=0, fix rates=0
/// @param target_CL Target lift coefficient
/// @return TrimConstraints configured for level flight
inline TrimConstraints LevelFlight(Scalar target_CL) {
    TrimConstraints constraints;

    // Alpha is free, will be adjusted to achieve CL
    constraints.alpha.type = ConstraintType::kCL;
    constraints.alpha.target = target_CL;

    // Beta fixed at zero
    constraints.beta.type = ConstraintType::kFixed;
    constraints.beta.target = 0.0;

    // All rates fixed at zero
    constraints.roll_rate.type = ConstraintType::kFixed;
    constraints.roll_rate.target = 0.0;
    constraints.pitch_rate.type = ConstraintType::kFixed;
    constraints.pitch_rate.target = 0.0;
    constraints.yaw_rate.type = ConstraintType::kFixed;
    constraints.yaw_rate.target = 0.0;

    return constraints;
}

/// Level flight with pitch trim: CL constraint + Cm=0
/// @param target_CL Target lift coefficient
/// @param elevator_name Name of elevator control surface
/// @return TrimConstraints configured for level trimmed flight
inline TrimConstraints LevelTrimmedFlight(Scalar target_CL,
                                          const std::string& elevator_name) {
    TrimConstraints constraints = LevelFlight(target_CL);

    // Add elevator to trim pitch moment to zero
    VariableConstraint elevator_constraint;
    elevator_constraint.type = ConstraintType::kPitchMoment;
    elevator_constraint.target = 0.0;
    constraints.control_constraints.emplace_back(elevator_name,
                                                 elevator_constraint);

    return constraints;
}

/// Steady turn: CL constraint + specified bank angle
/// @param target_CL Target lift coefficient
/// @param bank_angle Bank angle in radians (used for turn rate calculation)
/// @return TrimConstraints configured for steady turn
inline TrimConstraints SteadyTurn(Scalar target_CL,
                                  [[maybe_unused]] Scalar bank_angle) {
    TrimConstraints constraints;

    // Alpha adjusted to achieve CL
    constraints.alpha.type = ConstraintType::kCL;
    constraints.alpha.target = target_CL;

    // Beta fixed at zero
    constraints.beta.type = ConstraintType::kFixed;
    constraints.beta.target = 0.0;

    // Roll rate fixed at zero (steady turn, not rolling)
    constraints.roll_rate.type = ConstraintType::kFixed;
    constraints.roll_rate.target = 0.0;

    // Pitch rate will depend on turn
    constraints.pitch_rate.type = ConstraintType::kFree;

    // Yaw rate determined by turn
    constraints.yaw_rate.type = ConstraintType::kFree;

    return constraints;
}

/// Steady sideslip: fixed beta, CL constraint
/// @param target_CL Target lift coefficient
/// @param sideslip_angle Sideslip angle in radians
/// @return TrimConstraints configured for steady sideslip
inline TrimConstraints SteadySideslip(Scalar target_CL, Scalar sideslip_angle) {
    TrimConstraints constraints;

    // Alpha adjusted to achieve CL
    constraints.alpha.type = ConstraintType::kCL;
    constraints.alpha.target = target_CL;

    // Beta fixed at specified sideslip
    constraints.beta.type = ConstraintType::kFixed;
    constraints.beta.target = sideslip_angle;

    // All rates fixed at zero
    constraints.roll_rate.type = ConstraintType::kFixed;
    constraints.roll_rate.target = 0.0;
    constraints.pitch_rate.type = ConstraintType::kFixed;
    constraints.pitch_rate.target = 0.0;
    constraints.yaw_rate.type = ConstraintType::kFixed;
    constraints.yaw_rate.target = 0.0;

    return constraints;
}

}  // namespace TrimPresets

}  // namespace vlm

#endif  // VLM_TRIM_FLIGHT_CONSTRAINTS_H_
