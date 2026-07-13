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
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/viscous/coupler.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus::viscous {
namespace {

/// Analytic NACA 0012 coordinates in Selig order (upper TE → LE → lower
/// TE) with cosine spacing and a closed trailing edge.
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
        double x = 0.5 * (1.0 + std::cos(theta));
        coords(0, i) = x;
        coords(1, i) = thickness(x);
    }
    for (int i = 1; i < n_per_side; ++i) {
        double theta = kPi * i / (n_per_side - 1);
        double x = 0.5 * (1.0 - std::cos(theta));
        coords(0, n_per_side - 1 + i) = x;
        coords(1, n_per_side - 1 + i) = -thickness(x);
    }
    return coords;
}

Aircraft BuildSmallWing(Index n_spanwise) {
    const Scalar span = 10.0;
    const Scalar chord = 1.0;
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(static_cast<int>(n_spanwise));
    wing.SetChordwiseSpacing(SpacingType::kCosine);
    wing.SetSpanwiseSpacing(SpacingType::kCosine);
    Section root(Vec3(0.0, 0.0, 0.0), chord);
    Section tip(Vec3(0.0, span / 2.0, 0.0), chord);
    wing.AddSection(root);
    wing.AddSection(tip);
    Aircraft aircraft("ViscousCouplerBoraTest");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(span * chord);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);
    return aircraft;
}

std::vector<CouplingStation> UniformStations(const std::vector<Scalar>& etas,
                                             Scalar reynolds = 1.0e6) {
    Eigen::MatrixXd coords = Naca0012SeligCoords();
    std::vector<CouplingStation> stations;
    stations.reserve(etas.size());
    for (Scalar e : etas) {
        CouplingStation s;
        s.eta = e;
        s.coords = coords;
        s.reynolds = reynolds;
        stations.push_back(std::move(s));
    }
    return stations;
}

CouplerOptions BoraCouplerOptions(std::vector<CouplingStation> stations) {
    CouplerOptions opt;
    opt.stations = std::move(stations);
    opt.backend = AirfoilSolverBackend::kBora;
    opt.mach = 0.0;
    opt.max_iterations = 20;
    opt.tolerance = 1.0e-6;
    return opt;
}

TEST(ViscousCouplerBoraTest, ConvergesOnSmallWingAndAddsViscousDrag) {
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition flight;
    flight.alpha = 2.0 * kDegreesToRadians;

    ViscousVlmCoupler coupler(
        solver, flight, BoraCouplerOptions(UniformStations({0.05, 0.5, 0.95})));
    CouplerResult res = coupler.Solve();

    EXPECT_TRUE(res.converged) << "Final residual " << res.final_residual
                               << " after " << res.iterations << " iterations";
    EXPECT_GT(res.results.CDv, 0.0);
    EXPECT_LT(res.results.CDv, 0.05);
    EXPECT_GT(res.results.CD, res.results.CDi_nearfield);
    EXPECT_GT(res.results.CL, 0.05)
        << "wing must produce real lift, CL = " << res.results.CL;
}

TEST(ViscousCouplerBoraTest, CoupledGradientMatchesFullFiniteDifference) {
    auto BuildOptions = []() {
        CouplerOptions opt = BoraCouplerOptions(UniformStations({0.1, 0.9}));
        opt.max_iterations = 50;
        opt.tolerance = 1.0e-12;
        opt.max_step_rad = 5.0 * kDegreesToRadians;
        return opt;
    };

    const Scalar alpha_base = 2.0 * kDegreesToRadians;
    const Scalar h_alpha = 0.005 * kDegreesToRadians;

    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    FlightCondition flight;
    flight.alpha = alpha_base;

    ViscousVlmCoupler coupler(solver, flight, BuildOptions());
    CouplerResult base = coupler.Solve();
    ASSERT_TRUE(base.converged)
        << "base did not converge: iters=" << base.iterations
        << ", final_res=" << base.final_residual << ", CL=" << base.results.CL;

    CoupledGradients g = coupler.ComputeAlphaDerivatives();

    auto SolveAt = [&](Scalar alpha) -> CouplerResult {
        Aircraft a = BuildSmallWing(/*n_spanwise=*/3);
        VortexLatticeSolver s;
        s.SetAircraft(a);
        FlightCondition f;
        f.alpha = alpha;
        ViscousVlmCoupler c(s, f, BuildOptions());
        return c.Solve();
    };

    CouplerResult plus = SolveAt(alpha_base + h_alpha);
    CouplerResult minus = SolveAt(alpha_base - h_alpha);
    ASSERT_TRUE(plus.converged);
    ASSERT_TRUE(minus.converged);

    Scalar dCL_dalpha_fd =
        (plus.results.CL - minus.results.CL) / (2.0 * h_alpha);
    Scalar dCD_dalpha_fd =
        (plus.results.CD - minus.results.CD) / (2.0 * h_alpha);

    Scalar cl_rel_err = std::abs(g.dCL_dalpha - dCL_dalpha_fd) /
                        std::max(std::abs(dCL_dalpha_fd), Scalar(1e-3));
    Scalar cd_rel_err = std::abs(g.dCD_dalpha - dCD_dalpha_fd) /
                        std::max(std::abs(dCD_dalpha_fd), Scalar(1e-4));
    // Every link in the chain (VLM forward chain, network gradients) is
    // analytic, so agreement is limited by the FD reference's own
    // truncation/roundoff budget — the same 5e-6 floor as the Scirocco
    // coupled-gradient test.
    EXPECT_LT(cl_rel_err, 5.0e-6)
        << "analytic dCL/dα = " << g.dCL_dalpha
        << ", FD reference = " << dCL_dalpha_fd << ", rel err = " << cl_rel_err;
    EXPECT_LT(cd_rel_err, 5.0e-6)
        << "analytic dCD/dα = " << g.dCD_dalpha
        << ", FD reference = " << dCD_dalpha_fd << ", rel err = " << cd_rel_err;
}

TEST(ViscousCouplerBoraTest, ShapeDerivativesAreKulfanSized) {
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition flight;
    flight.alpha = 2.0 * kDegreesToRadians;

    std::vector<CouplingStation> stations = UniformStations({0.1, 0.9});
    ViscousVlmCoupler coupler(solver, flight, BoraCouplerOptions(stations));
    CouplerResult res = coupler.Solve();
    ASSERT_TRUE(res.converged);

    CoupledShapeGradients g = coupler.ComputeShapeDerivatives();
    ASSERT_EQ(g.dCL_dcoords.size(), stations.size());
    ASSERT_EQ(g.dCD_dcoords.size(), stations.size());
    for (std::size_t k = 0; k < stations.size(); ++k) {
        ASSERT_EQ(g.dCL_dcoords[k].rows(), kBoraNumShapeParams);
        ASSERT_EQ(g.dCL_dcoords[k].cols(), 1);
        ASSERT_EQ(g.dCD_dcoords[k].rows(), kBoraNumShapeParams);
        ASSERT_EQ(g.dCD_dcoords[k].cols(), 1);
        EXPECT_TRUE(g.dCL_dcoords[k].allFinite());
        EXPECT_TRUE(g.dCD_dcoords[k].allFinite());
    }
}

}  // namespace
}  // namespace pegasus::viscous
