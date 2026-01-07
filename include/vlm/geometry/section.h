// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_GEOMETRY_SECTION_H_
#define VLM_GEOMETRY_SECTION_H_

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "vlm/core/types.h"
#include "vlm/geometry/airfoil.h"

namespace vlm {

/// @brief Control surface definition on a wing section.
///
/// Defines a control surface (flap, aileron, elevator, etc.) that deflects
/// about a hinge line. The hinge is defined by chord fraction and hinge
/// vector direction.
struct ControlSurfaceDefinition {
  std::string name;        ///< Control surface name (e.g., "aileron")
  Scalar hinge_x_fraction; ///< Chordwise hinge location (0=LE, 1=TE)
  Vec3 hinge_vector;       ///< Hinge axis direction (typically (0,1,0) or (0,0,1))
  Scalar sign;             ///< Sign of deflection (+1 or -1)
  Scalar gain;             ///< Deflection gain multiplier
};

/// @brief Design variable definition for optimization.
///
/// Defines how a geometric parameter varies with a design variable,
/// used for sensitivity analysis and optimization.
struct DesignVariableDefinition {
  std::string name;  ///< Design variable name
  Scalar weight;     ///< Weighting factor for this section
};

/// @brief Wing section definition for a lifting surface.
///
/// A section defines a spanwise station on a lifting surface, specifying
/// the leading edge position, chord, twist angle, and airfoil. Multiple
/// sections define the planform of a surface through linear interpolation.
///
/// In AVL terminology, this corresponds to a "SECTION" block within a
/// "SURFACE" definition.
class Section {
 public:
  Section() = default;

  /// @brief Creates a section with basic parameters.
  ///
  /// @param leading_edge_position Position of the leading edge (x, y, z)
  /// @param chord Chord length
  /// @param twist_angle Angle of incidence in degrees (positive = LE up)
  /// @param airfoil Optional airfoil (defaults to flat plate)
  Section(const Vec3& leading_edge_position, Scalar chord,
          Scalar twist_angle = 0.0,
          std::shared_ptr<Airfoil> airfoil = nullptr);

  /// @name Basic Geometry
  /// @{

  /// @brief Gets the leading edge position.
  const Vec3& GetLeadingEdgePosition() const { return leading_edge_position_; }

  /// @brief Sets the leading edge position.
  void SetLeadingEdgePosition(const Vec3& position) {
    leading_edge_position_ = position;
  }

  /// @brief Gets the chord length.
  Scalar GetChord() const { return chord_; }

  /// @brief Sets the chord length.
  void SetChord(Scalar chord) { chord_ = chord; }

  /// @brief Gets the twist angle in degrees.
  Scalar GetTwistAngle() const { return twist_angle_deg_; }

  /// @brief Sets the twist angle in degrees.
  void SetTwistAngle(Scalar angle_deg) { twist_angle_deg_ = angle_deg; }

  /// @brief Gets the twist angle in radians.
  Scalar GetTwistAngleRadians() const {
    return twist_angle_deg_ * kDegreesToRadians;
  }

  /// @}

  /// @name Airfoil
  /// @{

  /// @brief Gets the airfoil.
  const Airfoil& GetAirfoil() const;

  /// @brief Gets the airfoil (non-const).
  Airfoil& GetAirfoilMutable();

  /// @brief Sets the airfoil.
  void SetAirfoil(std::shared_ptr<Airfoil> airfoil) { airfoil_ = airfoil; }

  /// @brief Checks if a custom airfoil is set.
  bool HasAirfoil() const { return airfoil_ != nullptr; }

  /// @}

  /// @name Aerodynamic Modifiers
  /// @{

  /// @brief Gets the CL-alpha scaling factor.
  ///
  /// Factor to scale the section lift curve slope. Default is 1.0.
  /// Values less than 1.0 can model flow separation or 3D effects.
  Scalar GetClAlphaFactor() const { return cl_alpha_factor_; }

  /// @brief Sets the CL-alpha scaling factor.
  void SetClAlphaFactor(Scalar factor) { cl_alpha_factor_ = factor; }

  /// @}

  /// @name Control Surfaces
  /// @{

  /// @brief Adds a control surface to this section.
  void AddControlSurface(const ControlSurfaceDefinition& control) {
    control_surfaces_.push_back(control);
  }

  /// @brief Gets all control surface definitions.
  const std::vector<ControlSurfaceDefinition>& GetControlSurfaces() const {
    return control_surfaces_;
  }

  /// @brief Gets the number of control surfaces.
  size_t GetNumControlSurfaces() const { return control_surfaces_.size(); }

  /// @}

  /// @name Design Variables
  /// @{

  /// @brief Adds a design variable to this section.
  void AddDesignVariable(const DesignVariableDefinition& design_var) {
    design_variables_.push_back(design_var);
  }

  /// @brief Gets all design variable definitions.
  const std::vector<DesignVariableDefinition>& GetDesignVariables() const {
    return design_variables_;
  }

  /// @}

  /// @name Paneling Overrides
  /// @{

  /// @brief Gets the spanwise panel count override.
  /// @return Panel count, or 0 if using surface default
  int GetSpanwisePanelCount() const { return spanwise_panel_count_; }

  /// @brief Sets the spanwise panel count override.
  void SetSpanwisePanelCount(int count) { spanwise_panel_count_ = count; }

  /// @brief Gets the spanwise spacing type.
  SpacingType GetSpanwiseSpacing() const { return spanwise_spacing_type_; }

  /// @brief Sets the spanwise spacing type.
  void SetSpanwiseSpacing(SpacingType spacing) { spanwise_spacing_type_ = spacing; }

  /// @brief Sets the spanwise spacing from a numeric parameter.
  void SetSpanwiseSpacing(Scalar spacing) {
    if (std::abs(spacing) < 0.01) {
      spanwise_spacing_type_ = SpacingType::kUniform;
    } else if (spacing > 0.99) {
      spanwise_spacing_type_ = SpacingType::kCosine;
    } else if (spacing < -0.99) {
      spanwise_spacing_type_ = SpacingType::kSine;
    } else if (spacing > 0) {
      spanwise_spacing_type_ = SpacingType::kCosine;
    } else {
      spanwise_spacing_type_ = SpacingType::kSine;
    }
  }

  /// @}

  /// @name Viscous Drag Model
  /// @{

  /// @brief Gets the viscous polar coefficients.
  ///
  /// 6-element array defining parabolic drag polar at this section.
  const std::array<Scalar, 6>& GetViscousPolar() const { return viscous_polar_; }

  /// @brief Sets the viscous polar coefficients.
  void SetViscousPolar(const std::array<Scalar, 6>& polar) { viscous_polar_ = polar; }

  /// @}

 private:
  Vec3 leading_edge_position_{0, 0, 0};  ///< Leading edge (x, y, z)
  Scalar chord_{1.0};                     ///< Chord length
  Scalar twist_angle_deg_{0.0};           ///< Angle of incidence (degrees)
  std::shared_ptr<Airfoil> airfoil_;      ///< Airfoil definition
  Scalar cl_alpha_factor_{1.0};           ///< CL-alpha scaling factor

  std::vector<ControlSurfaceDefinition> control_surfaces_;
  std::vector<DesignVariableDefinition> design_variables_;

  int spanwise_panel_count_{0};     ///< Override panel count (0 = use default)
  SpacingType spanwise_spacing_type_{SpacingType::kNone};  ///< Spacing type

  std::array<Scalar, 6> viscous_polar_{0, 0, 0, 0, 0, 0};  ///< Viscous polar
};

}  // namespace vlm

#endif  // VLM_GEOMETRY_SECTION_H_
