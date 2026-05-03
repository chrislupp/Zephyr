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

#ifndef PEGASUS_GEOMETRY_AIRFOIL_H_
#define PEGASUS_GEOMETRY_AIRFOIL_H_

#include "pegasus/core/types.h"
#include "pegasus/utils/cubic_spline.h"

#include <string>
#include <vector>

namespace pegasus {

/// @brief Represents an airfoil camber line for vortex lattice analysis.
///
/// In VLM, only the camber line (mean line) of the airfoil affects the
/// aerodynamic solution. Thickness is optionally stored for visualization
/// and viscous drag estimation.
///
/// Airfoils can be created from:
/// - NACA 4-digit designation (e.g., "2412")
/// - NACA 5-digit designation (e.g., "23012")
/// - Coordinate file (Selig or Lednicer format)
/// - Direct camber/thickness specification
class Airfoil {
public:
    Airfoil() = default;

    /// @brief Creates a flat plate airfoil (zero camber).
    /// @return Flat plate airfoil
    static Airfoil CreateFlatPlate();

    /// @brief Creates an airfoil from NACA 4-digit designation.
    ///
    /// NACA 4-digit format: MPTT
    /// - M: Maximum camber as percentage of chord (0-9)
    /// - P: Position of maximum camber in tenths of chord (0-9)
    /// - TT: Maximum thickness as percentage of chord (01-99)
    ///
    /// Example: "2412" = 2% camber at 40% chord, 12% thick
    ///
    /// @param naca_code 4-digit NACA designation string or integer
    /// @param num_points Number of points for camber line (default 50)
    /// @return Airfoil object
    static Airfoil CreateFromNaca4Digit(const std::string& naca_code,
                                        int num_points = 50);
    static Airfoil CreateFromNaca4Digit(int naca_code, int num_points = 50);

    /// @brief Creates an airfoil from coordinate file.
    ///
    /// Supports Selig format (coordinates starting from upper surface trailing
    /// edge, going around leading edge to lower surface trailing edge) and
    /// Lednicer format (upper and lower surfaces listed separately).
    ///
    /// @param filename Path to coordinate file
    /// @param num_camber_points Number of points for camber line interpolation
    /// @return Airfoil object
    static Airfoil CreateFromCoordinateFile(const std::string& filename,
                                            int num_camber_points = 50);

    /// @brief Creates an airfoil from explicit coordinate arrays.
    ///
    /// @param x_coords X coordinates (0 to 1, normalized by chord)
    /// @param y_coords Y coordinates (camber line or full airfoil)
    /// @param is_camber_line True if y_coords is camber line, false if full
    /// airfoil
    /// @return Airfoil object
    static Airfoil CreateFromCoordinates(const std::vector<Scalar>& x_coords,
                                         const std::vector<Scalar>& y_coords,
                                         bool is_camber_line = true);

    /// @brief Evaluates camber line height at a given chord location.
    /// @param x Chordwise location (0 = LE, 1 = TE)
    /// @return Camber line height (y/c)
    [[nodiscard]] Scalar EvaluateCamber(Scalar x) const;

    /// @brief Evaluates camber line slope (dy/dx) at a given chord location.
    /// @param x Chordwise location
    /// @return Camber line slope
    [[nodiscard]] Scalar EvaluateCamberSlope(Scalar x) const;

    /// @brief Evaluates thickness at a given chord location.
    /// @param x Chordwise location
    /// @return Thickness (t/c)
    [[nodiscard]] Scalar EvaluateThickness(Scalar x) const;

    /// @brief Returns the camber line x-coordinates.
    [[nodiscard]] const std::vector<Scalar>& GetCamberX() const {
        return camber_x_;
    }

    /// @brief Returns the camber line y-coordinates.
    [[nodiscard]] const std::vector<Scalar>& GetCamberY() const {
        return camber_y_;
    }

    /// @brief Returns the maximum camber value.
    [[nodiscard]] Scalar GetMaxCamber() const { return max_camber_; }

    /// @brief Returns the chordwise position of maximum camber.
    [[nodiscard]] Scalar GetMaxCamberPosition() const {
        return max_camber_position_;
    }

    /// @brief Returns the maximum thickness.
    [[nodiscard]] Scalar GetMaxThickness() const { return max_thickness_; }

    /// @brief Returns the airfoil name/designation.
    [[nodiscard]] const std::string& GetName() const { return name_; }

    /// @brief Sets the airfoil name.
    void SetName(const std::string& name) { name_ = name; }

    /// @brief Checks if the airfoil has been initialized.
    [[nodiscard]] bool IsValid() const { return !camber_x_.empty(); }

private:
    /// @brief Builds the camber spline from stored coordinates.
    void BuildCamberSpline();

    /// @brief Extracts camber and thickness from full airfoil coordinates.
    ///
    /// Given upper and lower surface coordinates, computes the mean camber
    /// line as (y_upper + y_lower) / 2 and thickness as (y_upper - y_lower).
    ///
    /// @param x_coords X coordinates of airfoil (may be Selig or Lednicer
    /// format)
    /// @param y_coords Y coordinates of airfoil
    void ExtractCamberFromCoordinates(const std::vector<Scalar>& x_coords,
                                      const std::vector<Scalar>& y_coords);

    std::string name_;                 ///< Airfoil designation
    std::vector<Scalar> camber_x_;     ///< Camber line x coordinates
    std::vector<Scalar> camber_y_;     ///< Camber line y coordinates
    std::vector<Scalar> thickness_;    ///< Thickness distribution
    CubicSpline camber_spline_;        ///< Spline for camber interpolation
    CubicSpline thickness_spline_;     ///< Spline for thickness interpolation
    Scalar max_camber_{0.0};           ///< Maximum camber value
    Scalar max_camber_position_{0.0};  ///< Position of max camber
    Scalar max_thickness_{0.0};        ///< Maximum thickness
};

}  // namespace pegasus

#endif  // PEGASUS_GEOMETRY_AIRFOIL_H_
