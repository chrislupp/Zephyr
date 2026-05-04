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

#pragma once

/// @file output_selection.h
/// @brief Selection of which aerodynamic outputs to compute.

#include <vector>

namespace pegasus::gradients {

/// @brief Selection of which aerodynamic outputs to compute
/// gradients for.
///
/// By default all outputs are enabled. Disabling unused outputs
/// reduces the number of adjoint solves and O(n^2)
/// velocity-kernel operations, yielding significant speedups
/// when only a subset is needed.
///
/// Example usage:
/// @code
/// OutputSelection sel;
/// sel.CY = false;
/// sel.Cl = false;
/// sel.Cm = false;
/// sel.Cn = false;
/// // Now only CL, CD, and CDi gradients are computed
/// calc.SetOutputSelection(sel);
/// @endcode
struct OutputSelection {
    bool CL = true;   ///< Lift coefficient
    bool CD = true;   ///< Drag coefficient
    bool CY = true;   ///< Side force coefficient
    bool Cl = true;   ///< Rolling moment coefficient
    bool Cm = true;   ///< Pitching moment coefficient
    bool Cn = true;   ///< Yawing moment coefficient
    bool CDi = true;  ///< Induced drag coefficient

    /// Number of enabled primary outputs (CL..Cn, excludes CDi).
    [[nodiscard]] int CountPrimary() const {
        return int(CL) + int(CD) + int(CY) + int(Cl) + int(Cm) + int(Cn);
    }

    /// Check if a primary output index (0-5) is enabled.
    [[nodiscard]] bool IsEnabled(int index) const {
        switch (index) {
            case 0:
                return CL;
            case 1:
                return CD;
            case 2:
                return CY;
            case 3:
                return Cl;
            case 4:
                return Cm;
            case 5:
                return Cn;
            default:
                return false;
        }
    }

    /// Build a list of enabled primary output indices.
    [[nodiscard]] std::vector<int> EnabledIndices() const {
        std::vector<int> indices;
        for (int i = 0; i < 6; ++i) {
            if (IsEnabled(i))
                indices.push_back(i);
        }
        return indices;
    }

    /// Return a selection with all outputs enabled.
    static OutputSelection All() { return OutputSelection{}; }
};

}  // namespace pegasus::gradients
