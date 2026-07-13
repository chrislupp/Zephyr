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

#ifndef PEGASUS_VISCOUS_AIRFOIL_SOLVER_H_
#define PEGASUS_VISCOUS_AIRFOIL_SOLVER_H_

#include "pegasus/core/types.h"

#include <Eigen/Core>

namespace pegasus::viscous {

/// @brief Selects the 2D viscous airfoil solver used at each coupling
/// station.
enum class AirfoilSolverBackend {
    /// Scirocco viscous panel + boundary-layer code (requires building with
    /// PEGASUS_ENABLE_SCIROCCO).
    kScirocco,
    /// Bora, a C++ port of the NeuralFoil surrogate of XFoil (requires
    /// building with PEGASUS_ENABLE_BORA).
    kBora,
};

/// @brief Abstract per-station 2D viscous airfoil solver.
///
/// The coupler owns one instance per coupling station and drives it in
/// fixed-α mode at every Newton iteration. All angles at this interface are
/// in radians; backends convert internally to their own conventions.
///
/// The shape-derivative layout (GetDclDx / GetDcdDx) is backend-defined:
/// Scirocco returns a 2×N matrix over its internal panel-node coordinates,
/// while Bora returns an 18×1 matrix over the Kulfan (CST) parameters
/// [8 upper weights, 8 lower weights, LE weight, TE thickness]. The coupled
/// adjoint in ViscousVlmCoupler::ComputeShapeDerivatives is agnostic to the
/// layout — it only scales and combines these matrices.
class AirfoilSolver {
public:
    virtual ~AirfoilSolver() = default;

    /// @brief Solve at fixed α (rad), Reynolds, Mach. Returns true on a
    /// converged (trustworthy) solution. On success, coefficient and
    /// sensitivity getters reflect the new state.
    virtual bool SolveFixedAlpha(Scalar alpha_rad, Scalar reynolds,
                                 Scalar mach) = 0;

    /// @brief Angle of attack of the most recent solve, radians.
    [[nodiscard]] virtual Scalar GetAlphaRadians() const = 0;

    /// @brief Convergence flag of the most recent solve.
    [[nodiscard]] virtual bool IsConverged() const = 0;

    /// @brief Section lift coefficient of the most recent solve.
    [[nodiscard]] virtual Scalar GetCl() const = 0;
    /// @brief Section drag coefficient of the most recent solve.
    [[nodiscard]] virtual Scalar GetCd() const = 0;
    /// @brief Section pitching moment coefficient of the most recent solve.
    [[nodiscard]] virtual Scalar GetCm() const = 0;

    /// @brief dCl/dα in units of per-radian.
    [[nodiscard]] virtual Scalar GetDclDalphaPerRadian() const = 0;
    /// @brief dCd/dα in units of per-radian.
    [[nodiscard]] virtual Scalar GetDcdDalphaPerRadian() const = 0;

    /// @brief dCl w.r.t. the backend's shape parameterization (see class
    /// docs for the layout).
    [[nodiscard]] virtual const Eigen::MatrixXd& GetDclDx() const = 0;
    /// @brief dCd w.r.t. the backend's shape parameterization — same layout
    /// as GetDclDx.
    [[nodiscard]] virtual const Eigen::MatrixXd& GetDcdDx() const = 0;
};

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_AIRFOIL_SOLVER_H_
