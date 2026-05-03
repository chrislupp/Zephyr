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

#ifndef PEGASUS_DYNAMICS_MASS_PROPERTIES_H_
#define PEGASUS_DYNAMICS_MASS_PROPERTIES_H_

#include "pegasus/core/types.h"

namespace pegasus {

/// Mass and inertia properties for dynamic analysis
///
/// Uses standard aircraft body-axis convention:
/// - X: forward (positive out the nose)
/// - Y: right wing (positive out right wing)
/// - Z: down (positive down)
///
/// Inertia tensor is defined as:
/// [Ixx  -Ixy  -Ixz]
/// [-Ixy  Iyy  -Iyz]
/// [-Ixz -Iyz   Izz]
///
/// Note: Products of inertia (Ixy, Iyz, Ixz) are stored as positive values
/// but appear with negative signs in the inertia tensor.
struct MassProperties {
    Scalar mass{0.0};  ///< Total mass

    Vec3 cg{Vec3::Zero()};  ///< Center of gravity position

    /// Moments of inertia about CG
    Scalar Ixx{0.0};  ///< Roll moment of inertia
    Scalar Iyy{0.0};  ///< Pitch moment of inertia
    Scalar Izz{0.0};  ///< Yaw moment of inertia

    /// Products of inertia about CG (stored as positive values)
    Scalar Ixy{0.0};  ///< XY product of inertia
    Scalar Iyz{0.0};  ///< YZ product of inertia
    Scalar Ixz{0.0};  ///< XZ product of inertia

    /// Get the full 3x3 inertia tensor
    [[nodiscard]] Mat3 GetInertiaTensor() const {
        Mat3 I;
        I << Ixx, -Ixy, -Ixz, -Ixy, Iyy, -Iyz, -Ixz, -Iyz, Izz;
        return I;
    }

    /// Set inertia from a 3x3 tensor
    void SetInertiaTensor(const Mat3& I) {
        Ixx = I(0, 0);
        Iyy = I(1, 1);
        Izz = I(2, 2);
        Ixy = -I(0, 1);
        Iyz = -I(1, 2);
        Ixz = -I(0, 2);
    }

    /// Check if mass properties are valid (positive mass and inertias)
    [[nodiscard]] bool IsValid() const {
        return mass > 0.0 && Ixx > 0.0 && Iyy > 0.0 && Izz > 0.0;
    }

    /// Get the inverse of the inertia tensor
    [[nodiscard]] Mat3 GetInverseInertiaTensor() const {
        return GetInertiaTensor().inverse();
    }
};

/// Flight state for dynamic analysis
///
/// Defines the complete state of the aircraft for eigenmode analysis.
/// Uses body-axis velocities and rates.
struct FlightState {
    Scalar velocity{0.0};  ///< Airspeed magnitude (m/s)
    Scalar alpha{0.0};     ///< Angle of attack (rad)
    Scalar beta{0.0};      ///< Sideslip angle (rad)

    /// Body-axis angular rates (non-dimensional: rate * length / V)
    Scalar roll_rate{0.0};   ///< Roll rate p*b/2V
    Scalar pitch_rate{0.0};  ///< Pitch rate q*c/2V
    Scalar yaw_rate{0.0};    ///< Yaw rate r*b/2V

    /// Euler angles (rad)
    Scalar phi{0.0};    ///< Bank angle
    Scalar theta{0.0};  ///< Pitch angle
    Scalar psi{0.0};    ///< Heading angle

    /// Atmospheric properties
    Scalar density{1.225};  ///< Air density (kg/m^3)
    Scalar gravity{9.81};   ///< Gravitational acceleration (m/s^2)

    /// Compute body-axis velocity components
    [[nodiscard]] Vec3 GetBodyVelocity() const {
        Scalar ca = std::cos(alpha);
        Scalar sa = std::sin(alpha);
        Scalar cb = std::cos(beta);
        Scalar sb = std::sin(beta);
        return velocity * Vec3(ca * cb, sb, sa * cb);
    }

    /// Compute dynamic pressure
    [[nodiscard]] Scalar GetDynamicPressure() const {
        return 0.5 * density * velocity * velocity;
    }
};

/// Apparent mass coefficients for added mass effects
///
/// Added mass arises from the inertia of air accelerated by the body.
/// For slender bodies and wings, this can be significant.
struct ApparentMass {
    /// 3x3 apparent mass matrix (adds to real mass for accelerations)
    Mat3 mass_matrix{Mat3::Zero()};

    /// 3x3 apparent inertia matrix (adds to real inertia)
    Mat3 inertia_matrix{Mat3::Zero()};

    /// Whether apparent mass effects are included
    bool include_effects{false};
};

}  // namespace pegasus

#endif  // PEGASUS_DYNAMICS_MASS_PROPERTIES_H_
