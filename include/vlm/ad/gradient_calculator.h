// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file gradient_calculator.h
/// @brief Gradient computation for aerodynamic optimization using Enzyme AD.
///
/// This module provides automatic differentiation of VLM aerodynamic outputs
/// with respect to flight conditions and geometry parameters, enabling
/// gradient-based optimization workflows.

#include <functional>

#include "vlm/ad/enzyme_interface.h"
#include "vlm/aero/aerodynamic_results.h"
#include "vlm/core/types.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/solver/vortex_lattice_solver.h"

namespace vlm {
namespace ad {

/// @brief Flight condition gradient structure.
///
/// Contains gradients of an objective function with respect to flight
/// condition parameters.
struct FlightConditionGradients {
  Scalar d_alpha{0.0};      ///< Gradient w.r.t. angle of attack
  Scalar d_beta{0.0};       ///< Gradient w.r.t. sideslip angle
  Scalar d_roll_rate{0.0};  ///< Gradient w.r.t. roll rate
  Scalar d_pitch_rate{0.0}; ///< Gradient w.r.t. pitch rate
  Scalar d_yaw_rate{0.0};   ///< Gradient w.r.t. yaw rate
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
  /// Gradients with respect to angle of attack (per radian)
  AeroGradients alpha_derivs;

  /// Gradients with respect to sideslip (per radian)
  AeroGradients beta_derivs;

  /// Gradients with respect to roll rate (non-dimensional)
  AeroGradients p_derivs;

  /// Gradients with respect to pitch rate (non-dimensional)
  AeroGradients q_derivs;

  /// Gradients with respect to yaw rate (non-dimensional)
  AeroGradients r_derivs;
};

// ============================================================================
// Geometry Gradient Structures
// ============================================================================

/// @brief Geometry design variable types for sensitivity analysis.
enum class GeometryVariable {
  kChord,    ///< Chord length at section
  kTwist,    ///< Twist angle at section (degrees)
  kXle,      ///< X-coordinate of leading edge
  kYle,      ///< Y-coordinate of leading edge (spanwise position)
  kZle       ///< Z-coordinate of leading edge (dihedral)
};

/// @brief Gradients of all aero outputs w.r.t. a geometry variable at one section.
///
/// These are computed using the adjoint method:
///   dJ/dθ = ∂J/∂θ|_Γ - λᵀ·(∂A/∂θ)·Γ + λᵀ·(∂b/∂θ)
///
/// where θ is a geometry parameter (chord, twist, position).
struct SectionGeometryGradients {
  std::string surface_name;  ///< Name of the surface containing this section
  Index section_index;       ///< Index of the section within the surface

  /// Gradients w.r.t. chord at this section
  AeroGradients chord_derivs;

  /// Gradients w.r.t. twist at this section (per degree)
  AeroGradients twist_derivs;

  /// Gradients w.r.t. x-position of leading edge
  AeroGradients xle_derivs;

  /// Gradients w.r.t. y-position (spanwise) of leading edge
  AeroGradients yle_derivs;

  /// Gradients w.r.t. z-position (vertical) of leading edge
  AeroGradients zle_derivs;
};

/// @brief Complete geometry gradients for an entire aircraft.
///
/// Contains gradients for all sections across all surfaces.
struct GeometryGradients {
  /// Per-section gradients, organized by surface then section
  std::vector<SectionGeometryGradients> section_gradients;

  /// Number of total design variables (5 per section: chord, twist, x, y, z)
  Index TotalDesignVariables() const { return section_gradients.size() * 5; }

  /// Get gradient for a specific coefficient w.r.t. all chord values
  VecX GetChordGradientsFor(const std::function<Scalar(const AeroGradients&)>& selector) const {
    VecX result(section_gradients.size());
    for (size_t i = 0; i < section_gradients.size(); ++i) {
      result(i) = selector(section_gradients[i].chord_derivs);
    }
    return result;
  }

  /// Get gradient for a specific coefficient w.r.t. all twist values
  VecX GetTwistGradientsFor(const std::function<Scalar(const AeroGradients&)>& selector) const {
    VecX result(section_gradients.size());
    for (size_t i = 0; i < section_gradients.size(); ++i) {
      result(i) = selector(section_gradients[i].twist_derivs);
    }
    return result;
  }
};

/// @brief Calculator for aerodynamic gradients using automatic differentiation.
///
/// This class provides methods to compute gradients of aerodynamic outputs
/// with respect to flight conditions. When compiled with Enzyme (VLM_USE_ENZYME),
/// gradients are computed using automatic differentiation. Otherwise, finite
/// differences are used as a fallback.
///
/// Example usage:
/// @code
/// GradientCalculator calc;
/// calc.SetAircraft(aircraft);
///
/// FlightCondition cond;
/// cond.alpha = 5.0 * kDegreesToRadians;
///
/// // Compute all flight condition derivatives
/// FullGradients grads = calc.ComputeFlightConditionDerivatives(cond);
///
/// // grads.alpha_derivs.dCL is dCL/dalpha
/// // grads.alpha_derivs.dCm is dCm/dalpha
/// @endcode
class GradientCalculator {
 public:
  GradientCalculator();
  ~GradientCalculator();

  /// @brief Set the aircraft configuration.
  ///
  /// This generates the mesh and prepares for analysis.
  /// @param aircraft Aircraft geometry
  void SetAircraft(const Aircraft& aircraft);

  /// @brief Set the reference point for moment calculations.
  void SetReferencePoint(const Vec3& ref);

  /// @brief Compute derivatives with respect to all flight conditions.
  ///
  /// Computes dCL/dalpha, dCm/dalpha, dCL/dbeta, etc. for all six
  /// flight condition parameters (alpha, beta, p, q, r).
  ///
  /// @param condition The operating point about which to linearize
  /// @return Structure containing all gradients
  FullGradients ComputeFlightConditionDerivatives(
      const FlightCondition& condition);

  /// @brief Compute gradient of a single objective with respect to alpha.
  ///
  /// Uses Enzyme reverse mode when available.
  /// @param condition Flight condition
  /// @return Gradient structure for alpha derivatives
  AeroGradients ComputeAlphaDerivatives(const FlightCondition& condition);

  /// @brief Compute gradient with respect to sideslip beta.
  AeroGradients ComputeBetaDerivatives(const FlightCondition& condition);

  /// @brief Compute gradient with respect to roll rate p.
  AeroGradients ComputeRollRateDerivatives(const FlightCondition& condition);

  /// @brief Compute gradient with respect to pitch rate q.
  AeroGradients ComputePitchRateDerivatives(const FlightCondition& condition);

  /// @brief Compute gradient with respect to yaw rate r.
  AeroGradients ComputeYawRateDerivatives(const FlightCondition& condition);

  // ==========================================================================
  // Geometry Gradients
  // ==========================================================================

  /// @brief Compute gradients with respect to all geometry parameters.
  ///
  /// For each section in the aircraft, computes gradients of all aerodynamic
  /// coefficients with respect to:
  /// - Chord length
  /// - Twist angle
  /// - Leading edge position (x, y, z)
  ///
  /// Uses the adjoint method:
  ///   dJ/dθ = ∂J/∂θ|_Γ - λᵀ·(∂A/∂θ)·Γ + λᵀ·(∂b/∂θ)
  ///
  /// @param condition Flight condition for analysis
  /// @return Structure containing gradients for all sections
  GeometryGradients ComputeGeometryGradients(const FlightCondition& condition);

  /// @brief Compute gradient w.r.t. a single geometry parameter at one section.
  ///
  /// @param condition Flight condition
  /// @param surface_index Index of the surface
  /// @param section_index Index of the section within the surface
  /// @param variable Which geometry variable (chord, twist, etc.)
  /// @return Aerodynamic coefficient gradients
  AeroGradients ComputeSectionGradient(
      const FlightCondition& condition,
      Index surface_index,
      Index section_index,
      GeometryVariable variable);

  /// @brief Check if using Enzyme AD (vs finite differences).
  static bool IsUsingEnzyme();

  /// @brief Set finite difference step size (for fallback mode).
  void SetFiniteDiffStep(Scalar step);

 private:
  /// Internal implementation details
  class Impl;
  std::unique_ptr<Impl> impl_;

  Scalar fd_step_{1e-6};
};

// ============================================================================
// Low-Level Enzyme Interface (Advanced Users)
// ============================================================================

#ifdef VLM_USE_ENZYME

/// @brief Pure function computing CL from alpha for Enzyme differentiation.
///
/// This function is designed to be differentiated by Enzyme. It takes
/// raw pointers to avoid STL containers that may confuse Enzyme.
///
/// @param alpha Angle of attack (radians)
/// @param circulation Pointer to circulation array [N]
/// @param control_points Pointer to control point coordinates [3N]
/// @param bound_vortex_left Pointer to left bound vortex coords [3N]
/// @param bound_vortex_right Pointer to right bound vortex coords [3N]
/// @param normals Pointer to panel normal vectors [3N]
/// @param strip_widths Pointer to strip widths [num_strips]
/// @param strip_chords Pointer to strip chords [num_strips]
/// @param strip_first_vortex Pointer to first vortex indices [num_strips]
/// @param strip_num_vortices Pointer to vortex counts per strip [num_strips]
/// @param num_vortices Total number of vortices
/// @param num_strips Total number of strips
/// @param reference_area Reference area for coefficient normalization
/// @param[out] result Pointer to output [CL, CD, Cm] array
void enzyme_compute_forces_kernel(
    double alpha,
    const double* circulation,
    const double* control_points,
    const double* bound_vortex_left,
    const double* bound_vortex_right,
    const double* normals,
    const double* strip_widths,
    const double* strip_chords,
    const int* strip_first_vortex,
    const int* strip_num_vortices,
    int num_vortices,
    int num_strips,
    double reference_area,
    double reference_chord,
    double* result);

/// @brief Differentiate force computation with respect to alpha.
///
/// @param alpha_grad Output: gradient dL/dalpha
/// @param circulation_grad Output: gradients dL/dGamma [N]
/// @param ... same as enzyme_compute_forces_kernel
void enzyme_compute_forces_gradient(
    double alpha, double* alpha_grad,
    const double* circulation, double* circulation_grad,
    const double* control_points,
    const double* bound_vortex_left,
    const double* bound_vortex_right,
    const double* normals,
    const double* strip_widths,
    const double* strip_chords,
    const int* strip_first_vortex,
    const int* strip_num_vortices,
    int num_vortices,
    int num_strips,
    double reference_area,
    double reference_chord,
    const double* result_grad);

#endif  // VLM_USE_ENZYME

}  // namespace ad
}  // namespace vlm
