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
#include "pegasus/viscous/bora_polar.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus::viscous {
namespace {

/// Analytic NACA 0012 coordinates in Selig order (upper TE → LE → lower
/// TE) with cosine spacing and a closed (zero-thickness) trailing edge.
/// Shape: 2×(2·n_per_side − 1).
Eigen::MatrixXd Naca0012SeligCoords(int n_per_side = 81) {
    auto thickness = [](double x) {
        return 5.0 * 0.12 *
               (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                0.2843 * x * x * x - 0.1036 * x * x * x * x);
    };
    const int n_points = 2 * n_per_side - 1;
    Eigen::MatrixXd coords(2, n_points);
    for (int i = 0; i < n_per_side; ++i) {
        double theta = kPi * i / (n_per_side - 1);
        double x = 0.5 * (1.0 + std::cos(theta));  // 1 → 0
        coords(0, i) = x;
        coords(1, i) = thickness(x);
    }
    for (int i = 1; i < n_per_side; ++i) {
        double theta = kPi * i / (n_per_side - 1);
        double x = 0.5 * (1.0 - std::cos(theta));  // 0 → 1
        coords(0, n_per_side - 1 + i) = x;
        coords(1, n_per_side - 1 + i) = -thickness(x);
    }
    return coords;
}

TEST(BoraPolarSolverTest, ProducesPlausibleNaca0012Coefficients) {
    BoraOptions options;
    BoraPolarSolver solver(Naca0012SeligCoords(), options);

    const Scalar alpha = 5.0 * kDegreesToRadians;
    ASSERT_TRUE(solver.SolveFixedAlpha(alpha, 1.0e6, 0.0));
    EXPECT_TRUE(solver.IsConverged());
    EXPECT_NEAR(solver.GetAlphaRadians(), alpha, 1e-15);

    // XFoil reference for NACA 0012 at α=5°, Re=1e6: cl ≈ 0.55, cd ≈ 0.007.
    EXPECT_GT(solver.GetCl(), 0.40);
    EXPECT_LT(solver.GetCl(), 0.75);
    EXPECT_GT(solver.GetCd(), 0.004);
    EXPECT_LT(solver.GetCd(), 0.02);
    // Symmetric section: quarter-chord moment stays small.
    EXPECT_LT(std::abs(solver.GetCm()), 0.05);
    // Lift slope in the linear range is close to 2π per radian.
    EXPECT_GT(solver.GetDclDalphaPerRadian(), 4.0);
    EXPECT_LT(solver.GetDclDalphaPerRadian(), 8.5);
    // The canonical airfoil must sit well inside the training distribution.
    EXPECT_GT(solver.GetAnalysisConfidence(), 0.9);
}

TEST(BoraPolarSolverTest, SymmetricSectionHasNearZeroLiftAtZeroAlpha) {
    BoraOptions options;
    BoraPolarSolver solver(Naca0012SeligCoords(), options);

    ASSERT_TRUE(solver.SolveFixedAlpha(0.0, 1.0e6, 0.0));
    EXPECT_LT(std::abs(solver.GetCl()), 0.02);
}

TEST(BoraPolarSolverTest, AcceptsClockwiseCoordinateOrder) {
    // Reversed columns = opposite loop orientation; the wrapper must detect
    // and normalize it, giving the same Kulfan fit and coefficients.
    const Eigen::MatrixXd selig = Naca0012SeligCoords();
    const Eigen::MatrixXd reversed = selig.rowwise().reverse();

    BoraOptions options;
    BoraPolarSolver a(selig, options);
    BoraPolarSolver b(reversed, options);

    const Scalar alpha = 3.0 * kDegreesToRadians;
    ASSERT_TRUE(a.SolveFixedAlpha(alpha, 1.0e6, 0.0));
    ASSERT_TRUE(b.SolveFixedAlpha(alpha, 1.0e6, 0.0));
    EXPECT_NEAR(a.GetCl(), b.GetCl(), 1e-12);
    EXPECT_NEAR(a.GetCd(), b.GetCd(), 1e-12);
}

TEST(BoraPolarSolverTest, AlphaGradientsMatchFiniteDifferences) {
    BoraOptions options;
    BoraPolarSolver solver(Naca0012SeligCoords(), options);

    const Scalar alpha = 4.0 * kDegreesToRadians;
    const Scalar re = 1.0e6;
    ASSERT_TRUE(solver.SolveFixedAlpha(alpha, re, 0.0));
    const Scalar dcl_dalpha = solver.GetDclDalphaPerRadian();
    const Scalar dcd_dalpha = solver.GetDcdDalphaPerRadian();

    const Scalar h = 0.01 * kDegreesToRadians;
    ASSERT_TRUE(solver.SolveFixedAlpha(alpha + h, re, 0.0));
    const Scalar cl_plus = solver.GetCl();
    const Scalar cd_plus = solver.GetCd();
    ASSERT_TRUE(solver.SolveFixedAlpha(alpha - h, re, 0.0));
    const Scalar cl_minus = solver.GetCl();
    const Scalar cd_minus = solver.GetCd();

    const Scalar dcl_fd = (cl_plus - cl_minus) / (2.0 * h);
    const Scalar dcd_fd = (cd_plus - cd_minus) / (2.0 * h);

    // The surrogate's gradients are analytic (forward-mode AD through the
    // network), so agreement is limited only by FD truncation error.
    EXPECT_NEAR(dcl_dalpha, dcl_fd,
                1e-4 * std::max(std::abs(dcl_fd), Scalar(1.0)));
    EXPECT_NEAR(dcd_dalpha, dcd_fd,
                1e-4 * std::max(std::abs(dcd_fd), Scalar(0.01)));
}

TEST(BoraPolarSolverTest, ShapeGradientsAreKulfanSizedAndFinite) {
    BoraOptions options;
    BoraPolarSolver solver(Naca0012SeligCoords(), options);
    ASSERT_TRUE(solver.SolveFixedAlpha(2.0 * kDegreesToRadians, 1.0e6, 0.0));

    const Eigen::MatrixXd& dcl_dx = solver.GetDclDx();
    const Eigen::MatrixXd& dcd_dx = solver.GetDcdDx();
    ASSERT_EQ(dcl_dx.rows(), kBoraNumShapeParams);
    ASSERT_EQ(dcl_dx.cols(), 1);
    ASSERT_EQ(dcd_dx.rows(), kBoraNumShapeParams);
    ASSERT_EQ(dcd_dx.cols(), 1);
    EXPECT_TRUE(dcl_dx.allFinite());
    EXPECT_TRUE(dcd_dx.allFinite());
    // Camber changes must move the lift: the gradient cannot vanish.
    EXPECT_GT(dcl_dx.cwiseAbs().maxCoeff(), 0.1);
}

TEST(BoraPolarSolverTest, SharedModelFactoryMatchesStandaloneSolver) {
    CouplingStation station;
    station.eta = 0.5;
    station.coords = Naca0012SeligCoords();
    station.reynolds = 1.0e6;

    BoraOptions options;
    auto shared = CreateBoraStationSolvers({station, station}, options);
    ASSERT_EQ(shared.size(), 2u);

    BoraPolarSolver standalone(station.coords, options);

    const Scalar alpha = 3.0 * kDegreesToRadians;
    ASSERT_TRUE(shared[0]->SolveFixedAlpha(alpha, station.reynolds, 0.0));
    ASSERT_TRUE(shared[1]->SolveFixedAlpha(alpha, station.reynolds, 0.0));
    ASSERT_TRUE(standalone.SolveFixedAlpha(alpha, station.reynolds, 0.0));

    EXPECT_NEAR(shared[0]->GetCl(), standalone.GetCl(), 1e-12);
    EXPECT_NEAR(shared[1]->GetCd(), standalone.GetCd(), 1e-12);
}

}  // namespace
}  // namespace pegasus::viscous
