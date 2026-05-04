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

#ifndef PEGASUS_SOLVER_HORSESHOE_VORTEX_H_
#define PEGASUS_SOLVER_HORSESHOE_VORTEX_H_

#include "pegasus/core/types.h"

namespace pegasus {

/// @brief Represents a horseshoe vortex element in the vortex lattice mesh.
///
/// A horseshoe vortex consists of:
/// - A bound vortex segment (from left to right endpoint)
/// - Two semi-infinite trailing vortices extending downstream
///
/// The vortex induces velocities according to the Biot-Savart law and
/// satisfies the flow tangency condition at its control point.
///
/// Reference: Katz & Plotkin, "Low-Speed Aerodynamics", Chapter 12
struct HorseshoeVortex {
    /// Left endpoint of bound vortex (edge 1)
    Vec3 bound_vortex_left{0, 0, 0};

    /// Right endpoint of bound vortex (edge 2)
    Vec3 bound_vortex_right{0, 0, 0};

    /// Center of bound vortex (midpoint between endpoints)
    Vec3 bound_vortex_center{0, 0, 0};

    /// Control point where flow tangency is enforced
    Vec3 control_point{0, 0, 0};

    /// Unit normal vector at control point
    Vec3 normal_at_control_point{0, 0, 1};

    /// Unit normal vector at bound vortex center
    Vec3 normal_at_vortex{0, 0, 1};

    /// Camber slope (dz/dx) at control point
    Scalar camber_slope_at_control_point{0.0};

    /// Camber slope (dz/dx) at bound vortex
    Scalar camber_slope_at_vortex{0.0};

    /// Chordwise width of this element (delta x/c * chord)
    Scalar chordwise_width{0.0};

    /// Chord length of the strip containing this element
    Scalar local_chord{0.0};

    /// Index of the strip containing this element
    Index strip_index{0};

    /// Index of the surface containing this element
    Index surface_index{0};

    /// Index of the logical component (for grouped surfaces)
    Index component_index{0};

    /// Whether this element responds to freestream alpha/beta
    bool affected_by_freestream{true};

    /// Whether flow tangency (V·n = 0) is enforced at this control point
    bool enforce_flow_tangency{true};

    /// Circulation strength (solved quantity)
    Scalar circulation{0.0};

    /// Pressure coefficient difference across element
    Scalar delta_cp{0.0};
};

/// @brief Represents a spanwise strip in the vortex lattice mesh.
///
/// A strip spans from one chordwise row of panels to the next and contains
/// multiple horseshoe vortex elements along the chord direction.
struct Strip {
    /// Leading edge position at left (inner) edge of strip
    Vec3 leading_edge_left{0, 0, 0};

    /// Leading edge position at right (outer) edge of strip
    Vec3 leading_edge_right{0, 0, 0};

    /// Leading edge position at strip center
    Vec3 leading_edge_center{0, 0, 0};

    /// Chord length at left edge
    Scalar chord_left{0.0};

    /// Chord length at right edge
    Scalar chord_right{0.0};

    /// Chord length at center
    Scalar chord_center{0.0};

    /// Strip width in y-z plane
    Scalar width{0.0};

    /// Leading edge sweep slope (dx_LE / dy)
    Scalar leading_edge_sweep_slope{0.0};

    /// Trailing edge sweep slope (dx_TE / dy)
    Scalar trailing_edge_sweep_slope{0.0};

    /// Local incidence (twist) angle at strip center [radians]
    Scalar twist_angle{0.0};

    /// Local incidence angle at left edge [radians]
    Scalar twist_angle_left{0.0};

    /// Local incidence angle at right edge [radians]
    Scalar twist_angle_right{0.0};

    /// Y-component of strip normal in Trefftz plane (for wake analysis)
    Scalar trefftz_normal_y{0.0};

    /// Z-component of strip normal in Trefftz plane (for wake analysis)
    Scalar trefftz_normal_z{1.0};

    /// Spanwise unit vector for velocity perpendicular definition
    Vec3 spanwise_unit_vector{0, 1, 0};

    /// Reference point for strip (used in force calculations)
    Vec3 reference_point{0, 0, 0};

    /// Index of first horseshoe vortex element in this strip
    Index first_vortex_index{0};

    /// Number of horseshoe vortex elements in this strip
    Index num_vortices{0};

    /// Index of the surface containing this strip
    Index surface_index{0};

    /// Viscous polar parameters [CL0, CD0, CL1, CD1, CL2, CD2]
    /// Used for profile drag estimation: CD = CD0 + dCD_dCL^2 * (CL - CL0)^2
    std::array<Scalar, 6> viscous_polar{0, 0, 0, 0, 0, 0};

    /// Whether viscous drag data is available for this strip
    bool has_viscous_data{false};

    /// Whether this strip is active (inside the fluid domain)
    bool is_active{true};
};

/// @brief Control surface influence on a vortex element.
struct ControlSurfaceInfluence {
    /// Name of the control surface
    std::string name;

    /// Gain factor (0-1) indicating fraction of element on control surface
    /// Combined with deflection to get effective normal rotation
    Scalar gain{0.0};

    /// Hinge axis unit vector
    Vec3 hinge_axis{0, 1, 0};

    /// Point on the hinge line (for moment calculations)
    Vec3 hinge_point{0, 0, 0};

    /// Reflection sign for y-duplicated surfaces (+1 or -1)
    Scalar reflection_sign{1.0};
};

}  // namespace pegasus

#endif  // PEGASUS_SOLVER_HORSESHOE_VORTEX_H_
