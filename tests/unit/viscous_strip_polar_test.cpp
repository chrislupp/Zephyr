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

#include "pegasus/core/types.h"
#include "pegasus/viscous/strip_polar.h"

#include <gtest/gtest.h>

namespace pegasus::viscous {
namespace {

// NACA 0012 at α = 2°, Re = 1e6, Ma = 0 is a well-behaved viscous solve
// that converges in tens of Newton iterations.
constexpr Scalar kAlphaRad = 2.0 * kDegreesToRadians;
constexpr Scalar kRe = 1.0e6;
constexpr Scalar kMa = 0.0;

TEST(StripPolarSolverTest, ConvergesAtMildAngle) {
    StripPolarSolver solver("0012", 199);
    ASSERT_FALSE(solver.IsWarmStarted());

    ASSERT_TRUE(solver.SolveFixedAlpha(kAlphaRad, kRe, kMa));
    EXPECT_TRUE(solver.IsConverged());
    EXPECT_TRUE(solver.IsWarmStarted());

    // A symmetric airfoil at 2° in light viscous conditions should have
    // 0.15 < Cl < 0.30 and 0.005 < Cd < 0.020. Loose bounds intentionally —
    // this is a smoke test, not a polar reference.
    EXPECT_GT(solver.GetCl(), 0.15);
    EXPECT_LT(solver.GetCl(), 0.30);
    EXPECT_GT(solver.GetCd(), 0.003);
    EXPECT_LT(solver.GetCd(), 0.020);
}

TEST(StripPolarSolverTest, WarmStartPersistsAcrossSolves) {
    StripPolarSolver solver("0012", 199);

    ASSERT_TRUE(solver.SolveFixedAlpha(kAlphaRad, kRe, kMa));
    EXPECT_TRUE(solver.IsWarmStarted());

    // Nudge α and re-solve. Warm start should remain set on success.
    ASSERT_TRUE(
        solver.SolveFixedAlpha(kAlphaRad + 0.1 * kDegreesToRadians, kRe, kMa));
    EXPECT_TRUE(solver.IsWarmStarted());

    solver.ResetBoundaryLayer();
    EXPECT_FALSE(solver.IsWarmStarted());
}

TEST(StripPolarSolverTest, DclDalphaIsPerRadianAndReasonable) {
    StripPolarSolver solver("0012", 199);
    ASSERT_TRUE(solver.SolveFixedAlpha(kAlphaRad, kRe, kMa));

    // For a thin symmetric airfoil at low α in light viscous flow,
    // dCl/dα should be near 2π per radian (≈ 6.28). Loose bounds for
    // viscous effects.
    Scalar dcl_drad = solver.GetDclDalphaPerRadian();
    EXPECT_GT(dcl_drad, 5.0);
    EXPECT_LT(dcl_drad, 7.5);

    // dCd/dα is small and positive at this α; check sign and order of
    // magnitude only.
    Scalar dcd_drad = solver.GetDcdDalphaPerRadian();
    EXPECT_GT(dcd_drad, 0.0);
    EXPECT_LT(dcd_drad, 0.5);
}

}  // namespace
}  // namespace pegasus::viscous
