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

#ifndef PEGASUS_VISCOUS_STRIP_POLAR_H_
#define PEGASUS_VISCOUS_STRIP_POLAR_H_

#include "pegasus/core/types.h"

#include <memory>
#include <string>

#include <Eigen/Core>

namespace scirocco {
class Scirocco;
}

namespace pegasus::viscous {

/// @brief Per-strip 2D viscous airfoil solver wrapping a Scirocco instance.
///
/// Owns a single Scirocco analysis object for one spanwise strip and
/// maintains warm-start state across coupling Newton iterations: the first
/// SolveFixedAlpha() initializes the boundary layer; subsequent solves
/// reuse the previously converged BL state by setting `oper.initbl = false`
/// on Scirocco.
///
/// All α inputs/outputs are in radians at this interface; the wrapper
/// internally converts to/from Scirocco's degree convention.
class StripPolarSolver {
public:
    /// @brief Construct from a NACA 4-digit designator.
    explicit StripPolarSolver(const std::string& naca = "0012",
                              int num_panels = 199);

    /// @brief Construct from a 2×N coordinate matrix (Scirocco order).
    ///
    /// `coords` is a 2-by-N matrix of airfoil node coordinates passed
    /// directly to Scirocco. `num_panels` controls Scirocco's internal
    /// repaneling.
    StripPolarSolver(const Eigen::MatrixXd& coords, int num_panels);

    ~StripPolarSolver();
    StripPolarSolver(StripPolarSolver&&) noexcept;
    StripPolarSolver& operator=(StripPolarSolver&&) noexcept;
    StripPolarSolver(const StripPolarSolver&) = delete;
    StripPolarSolver& operator=(const StripPolarSolver&) = delete;

    /// @brief Solve at fixed α (rad), Reynolds, Mach. Returns true on
    /// converged solution.
    ///
    /// On success, sensitivities are pre-computed and the warm-start state
    /// is cached for the next call. On failure, the warm-start state is
    /// dropped so that the next call performs a fresh BL initialization.
    bool SolveFixedAlpha(Scalar alpha_rad, Scalar reynolds, Scalar mach);

    /// @brief Convergence flag for the most recent solve.
    [[nodiscard]] bool IsConverged() const;

    /// @brief Whether the next solve will reuse a cached BL state.
    [[nodiscard]] bool IsWarmStarted() const { return warm_started_; }

    /// @brief Drop the cached BL state so the next solve does a fresh
    /// initialization.
    void ResetBoundaryLayer() { warm_started_ = false; }

    /// @brief Section lift coefficient from the most recent converged solve.
    [[nodiscard]] Scalar GetCl() const;
    /// @brief Section drag coefficient from the most recent converged solve.
    [[nodiscard]] Scalar GetCd() const;
    /// @brief Section pitching moment coefficient from the most recent
    /// converged solve.
    [[nodiscard]] Scalar GetCm() const;

    /// @brief dCl/dα in units of per-radian.
    [[nodiscard]] Scalar GetDclDalphaPerRadian() const;
    /// @brief dCd/dα in units of per-radian.
    [[nodiscard]] Scalar GetDcdDalphaPerRadian() const;

    /// @brief dCl / d(panel-node-coords). Returns a 2×N matrix matching the
    /// Scirocco-internal panel layout (which may differ from input coords if
    /// repaneling occurred).
    [[nodiscard]] const Eigen::MatrixXd& GetDclDx() const;
    /// @brief dCd / d(panel-node-coords) — same layout as GetDclDx.
    [[nodiscard]] const Eigen::MatrixXd& GetDcdDx() const;

    /// @brief Number of panels in the Scirocco-internal airfoil. May differ
    /// from `num_panels` passed at construction if Scirocco repaneled.
    [[nodiscard]] int GetNumPanels() const;

private:
    std::unique_ptr<scirocco::Scirocco> solver_;
    bool warm_started_{false};
};

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_STRIP_POLAR_H_
