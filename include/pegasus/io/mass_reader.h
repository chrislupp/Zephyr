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

#ifndef PEGASUS_IO_MASS_READER_H_
#define PEGASUS_IO_MASS_READER_H_

#include "pegasus/core/types.h"
#include "pegasus/dynamics/mass_properties.h"

#include <string>
#include <vector>

namespace pegasus {

/// @brief A single mass component from a .mass file.
struct MassComponent {
    std::string name;             ///< Component name/description
    Scalar mass{0.0};             ///< Mass of component
    Vec3 position{Vec3::Zero()};  ///< CG position (x, y, z)
    Scalar Ixx{0.0};              ///< Moment of inertia about x-axis
    Scalar Iyy{0.0};              ///< Moment of inertia about y-axis
    Scalar Izz{0.0};              ///< Moment of inertia about z-axis
    Scalar Ixy{0.0};              ///< Product of inertia
    Scalar Ixz{0.0};              ///< Product of inertia
    Scalar Iyz{0.0};              ///< Product of inertia
};

/// @brief Result from parsing a .mass file.
struct MassFileData {
    /// Unit conversion factors
    Scalar length_unit{
        1.0};  ///< Length unit scale (e.g., 0.0254 for inches to meters)
    Scalar mass_unit{1.0};  ///< Mass unit scale (e.g., 0.001 for grams to kg)
    Scalar time_unit{1.0};  ///< Time unit scale

    /// Default environment values
    Scalar gravity{9.81};   ///< Gravitational acceleration
    Scalar density{1.225};  ///< Air density

    /// Individual mass components
    std::vector<MassComponent> components;

    /// @brief Computes total mass properties from all components.
    ///
    /// Combines all components using parallel axis theorem to compute
    /// total mass, CG location, and inertia tensor about the CG.
    ///
    /// @return Combined mass properties in SI units
    [[nodiscard]] MassProperties ComputeTotalProperties() const;

    /// @brief Gets the total mass (sum of all components).
    [[nodiscard]] Scalar GetTotalMass() const;

    /// @brief Gets the center of gravity position.
    [[nodiscard]] Vec3 GetCenterOfGravity() const;
};

/// @brief Parser for AVL-format mass property files (.mass).
///
/// The mass file format defines:
/// - Unit conversions (Lunit, Munit, Tunit)
/// - Default gravity and density values
/// - Individual mass components with position and inertia
///
/// Example file format:
/// @code
/// Lunit = 1.0 m
/// Munit = 1.0 kg
/// Tunit = 1.0 s
/// g = 9.81
/// rho = 1.225
/// #  mass    x     y     z      Ixx    Iyy   Izz   Ixy  Ixz  Iyz
///   100.0   0.5   0.0   0.0    10.0   50.0  55.0    0    0    0  ! fuselage
///    50.0   0.3   2.0   0.1     5.0    1.0   5.0    0    0    0  ! right wing
/// @endcode
class MassReader {
public:
    MassReader() = default;

    /// @brief Reads mass properties from a file.
    ///
    /// @param filename Path to the .mass file
    /// @return Parsed mass file data
    /// @throws std::runtime_error if file cannot be opened or parsed
    [[nodiscard]] MassFileData ReadFile(const std::string& filename) const;

    /// @brief Reads mass properties from a string.
    ///
    /// @param content The mass file content as a string
    /// @return Parsed mass file data
    /// @throws std::runtime_error if parsing fails
    [[nodiscard]] MassFileData ReadString(const std::string& content) const;

private:
    /// Parses a single line from the mass file.
    void ParseLine(const std::string& line, MassFileData& data,
                   Vec3& scale_factors, Vec3& offset_factors) const;

    /// Parses a unit definition line (Lunit, Munit, Tunit).
    bool ParseUnitLine(const std::string& line, MassFileData& data) const;

    /// Parses a constant definition line (g, rho).
    bool ParseConstantLine(const std::string& line, MassFileData& data) const;

    /// Parses a scale factor line (starting with *).
    bool ParseScaleLine(const std::string& line, Vec3& scale_factors) const;

    /// Parses an offset line (starting with +).
    bool ParseOffsetLine(const std::string& line, Vec3& offset_factors) const;

    /// Parses a mass component line.
    bool ParseComponentLine(const std::string& line, MassFileData& data,
                            const Vec3& scale_factors,
                            const Vec3& offset_factors) const;
};

}  // namespace pegasus

#endif  // PEGASUS_IO_MASS_READER_H_
