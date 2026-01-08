// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_AERO_AERODYNAMIC_RESULTS_H_
#define VLM_AERO_AERODYNAMIC_RESULTS_H_

#include "vlm/core/types.h"

namespace vlm {

/// Force and moment coefficients for a single strip (spanwise section)
struct StripLoads {
  Index strip_index{0};

  /// Strip geometry
  Scalar chord{0.0};
  Scalar width{0.0};
  Scalar area{0.0};
  Vec3 leading_edge_position{Vec3::Zero()};
  Vec3 quarter_chord_position{Vec3::Zero()};

  /// Force coefficients (normalized by strip area)
  Scalar cl{0.0};   ///< Lift coefficient
  Scalar cd{0.0};   ///< Drag coefficient
  Scalar cy{0.0};   ///< Side force coefficient
  Vec3 cf{Vec3::Zero()};  ///< Force coefficients (CFx, CFy, CFz)

  /// Moment coefficients about strip quarter chord (normalized by strip chord*area)
  Vec3 cm{Vec3::Zero()};  ///< Moment coefficients (CMx, CMy, CMz)

  /// Span loading (c * Cn)
  Scalar cnc{0.0};

  /// Induced angle and velocity at control point
  Scalar alpha_induced{0.0};
  Vec3 induced_velocity{Vec3::Zero()};
};

/// Force and moment coefficients for a lifting surface
struct SurfaceLoads {
  Index surface_index{0};
  std::string name;

  /// Force coefficients (normalized by reference area)
  Scalar CL{0.0};   ///< Lift coefficient
  Scalar CD{0.0};   ///< Drag coefficient (near-field)
  Scalar CY{0.0};   ///< Side force coefficient
  Vec3 CF{Vec3::Zero()};  ///< Force coefficients

  /// Moment coefficients about reference point (normalized by ref area * ref length)
  Vec3 CM{Vec3::Zero()};  ///< Moment coefficients (roll, pitch, yaw)

  /// Surface area
  Scalar area{0.0};

  /// Strip-wise loads for this surface
  std::vector<StripLoads> strip_loads;
};

/// Complete aerodynamic results for an aircraft configuration
struct AerodynamicResults {
  bool is_valid{false};

  // =========================================================================
  // Reference quantities
  // =========================================================================
  Scalar reference_area{1.0};
  Scalar reference_chord{1.0};
  Scalar reference_span{1.0};
  Vec3 moment_reference_point{Vec3::Zero()};

  // =========================================================================
  // Total force coefficients (body axes)
  // =========================================================================

  /// Total body-axis force coefficients
  Scalar CX{0.0};   ///< Axial force (positive forward)
  Scalar CY{0.0};   ///< Side force (positive right)
  Scalar CZ{0.0};   ///< Normal force (positive down)
  Vec3 CF{Vec3::Zero()};  ///< Force vector (CX, CY, CZ)

  /// Total body-axis moment coefficients about reference point
  Scalar Cl{0.0};   ///< Rolling moment (positive right wing down)
  Scalar Cm{0.0};   ///< Pitching moment (positive nose up)
  Scalar Cn{0.0};   ///< Yawing moment (positive nose right)
  Vec3 CM{Vec3::Zero()};  ///< Moment vector (Cl, Cm, Cn)

  // =========================================================================
  // Stability axis coefficients
  // =========================================================================

  /// Lift coefficient (force perpendicular to freestream)
  Scalar CL{0.0};

  /// Drag coefficient (force along freestream, near-field)
  Scalar CD{0.0};

  /// Induced drag from near-field integration
  Scalar CDi_nearfield{0.0};

  /// Viscous/parasitic drag (from profile drag input)
  Scalar CDv{0.0};

  // =========================================================================
  // Far-field (Trefftz plane) results
  // =========================================================================

  /// Far-field lift coefficient
  Scalar CL_ff{0.0};

  /// Far-field side force coefficient
  Scalar CY_ff{0.0};

  /// Far-field induced drag coefficient
  Scalar CDi_ff{0.0};

  /// Span efficiency factor (Oswald efficiency)
  /// e = CL^2 / (pi * AR * CDi)
  Scalar span_efficiency{0.0};

  // =========================================================================
  // Surface and strip breakdown
  // =========================================================================

  std::vector<SurfaceLoads> surface_loads;

  // =========================================================================
  // Helper methods
  // =========================================================================

  /// Get aspect ratio from reference dimensions
  Scalar GetAspectRatio() const {
    return reference_span * reference_span / reference_area;
  }

  /// Get induced drag coefficient from ideal elliptic loading
  Scalar GetIdealInducedDrag() const {
    Scalar AR = GetAspectRatio();
    return CL * CL / (kPi * AR);
  }

  /// Reset all values to zero
  void Clear() {
    is_valid = false;
    CX = CY = CZ = 0.0;
    CF = Vec3::Zero();
    Cl = Cm = Cn = 0.0;
    CM = Vec3::Zero();
    CL = CD = CDi_nearfield = CDv = 0.0;
    CL_ff = CY_ff = CDi_ff = 0.0;
    span_efficiency = 0.0;
    surface_loads.clear();
  }
};

}  // namespace vlm

#endif  // VLM_AERO_AERODYNAMIC_RESULTS_H_
