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

#include <gtest/gtest.h>
#include <scirocco/scirocco.h>

namespace pegasus::viscous {
namespace {

// Pull a canonical NACA 0012 panel-coordinate matrix out of Scirocco so the
// station-based tests can reuse the same airfoil at multiple stations.
Eigen::MatrixXd Naca0012Coords(int num_panels = 99) {
    scirocco::Scirocco airfoil("0012", num_panels);
    return airfoil.GetFoilCoords();
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
    Aircraft aircraft("ViscousCouplerTest");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(span * chord);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);
    return aircraft;
}

/// Half-wing with leading-edge sweep: tip LE shifted aft by `sweep_x`.
Aircraft BuildSweptWing(Index n_spanwise, Scalar sweep_x) {
    const Scalar span = 10.0;
    const Scalar chord = 1.0;
    Surface wing("Wing");
    wing.SetChordwisePanelCount(4);
    wing.SetSpanwisePanelCount(static_cast<int>(n_spanwise));
    wing.SetChordwiseSpacing(SpacingType::kCosine);
    wing.SetSpanwiseSpacing(SpacingType::kCosine);
    Section root(Vec3(0.0, 0.0, 0.0), chord);
    Section tip(Vec3(sweep_x, span / 2.0, 0.0), chord);
    wing.AddSection(root);
    wing.AddSection(tip);
    Aircraft aircraft("ViscousCouplerSweptTest");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(span * chord);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);
    return aircraft;
}

// Build a list of coupling stations at user-specified eta values, all with
// the same NACA 0012 airfoil and Reynolds. Models the "uniform airfoil"
// simple case while exercising the multi-station path.
std::vector<CouplingStation> UniformStations(const std::vector<Scalar>& etas,
                                             int panels = 99,
                                             Scalar reynolds = 1.0e6) {
    Eigen::MatrixXd coords = Naca0012Coords(panels);
    std::vector<CouplingStation> stations;
    stations.reserve(etas.size());
    for (Scalar e : etas) {
        CouplingStation s;
        s.eta = e;
        s.coords = coords;
        s.num_panels = panels;
        s.reynolds = reynolds;
        stations.push_back(std::move(s));
    }
    return stations;
}

TEST(ViscousCouplerTest, ConvergesOnSmallWingAndAddsViscousDrag) {
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition flight;
    flight.alpha = 2.0 * kDegreesToRadians;

    CouplerOptions opt;
    opt.stations = UniformStations({0.05, 0.5, 0.95}, /*panels=*/99);
    opt.mach = 0.0;
    opt.max_iterations = 20;
    opt.tolerance = 1.0e-6;

    ViscousVlmCoupler coupler(solver, flight, opt);
    CouplerResult res = coupler.Solve();

    EXPECT_TRUE(res.converged) << "Final residual " << res.final_residual
                               << " after " << res.iterations << " iterations";
    EXPECT_GT(res.results.CDv, 0.0);
    EXPECT_LT(res.results.CDv, 0.05);
    EXPECT_GT(res.results.CD, res.results.CDi_nearfield);
    EXPECT_GT(res.results.CL, 0.05)
        << "wing must produce real lift, CL = " << res.results.CL;
}

TEST(ViscousCouplerTest, CoupledGradientMatchesFullFiniteDifference) {
    auto BuildOptions = []() {
        CouplerOptions opt;
        opt.stations = UniformStations({0.1, 0.9}, /*panels=*/99);
        opt.mach = 0.0;
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
    // Analytical chain agrees with the fully-re-converged FD reference to
    // ~3e-6 relative; the floor is set by the FD reference's intrinsic
    // central-difference roundoff/truncation budget for typical h, not by
    // any limit of the analytical chain. 5e-6 leaves headroom.
    EXPECT_LT(cl_rel_err, 5.0e-6)
        << "analytic dCL/dα = " << g.dCL_dalpha
        << ", FD reference = " << dCL_dalpha_fd << ", rel err = " << cl_rel_err;
    EXPECT_LT(cd_rel_err, 5.0e-6)
        << "analytic dCD/dα = " << g.dCD_dalpha
        << ", FD reference = " << dCD_dalpha_fd << ", rel err = " << cd_rel_err;
}

// Validate the coupled adjoint for per-station airfoil-coordinate DVs by
// FD-perturbing one coordinate, re-converging the coupler from scratch,
// and central-differencing CL and CD. Picks a single (station, coord)
// pair to keep the test fast — fully validating all 2·N entries per
// station would take hundreds of full coupler solves.
TEST(ViscousCouplerTest, ShapeGradientMatchesFiniteDifference) {
    auto BuildBaseStations = []() {
        return UniformStations({0.1, 0.9}, /*panels=*/99);
    };
    auto BuildOptions = [&](const std::vector<CouplingStation>& stations) {
        CouplerOptions opt;
        opt.stations = stations;
        opt.mach = 0.0;
        opt.max_iterations = 50;
        opt.tolerance = 1.0e-12;
        opt.max_step_rad = 5.0 * kDegreesToRadians;
        return opt;
    };

    const Scalar alpha_base = 2.0 * kDegreesToRadians;
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);

    std::vector<CouplingStation> base_stations = BuildBaseStations();
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    FlightCondition flight;
    flight.alpha = alpha_base;

    ViscousVlmCoupler coupler(solver, flight, BuildOptions(base_stations));
    CouplerResult base = coupler.Solve();
    ASSERT_TRUE(base.converged);

    CoupledShapeGradients g = coupler.ComputeShapeDerivatives();
    ASSERT_EQ(g.dCL_dcoords.size(), base_stations.size());
    ASSERT_EQ(g.dCD_dcoords.size(), base_stations.size());

    // Pick station 0, the z-coordinate of a mid-upper-surface panel node.
    // Upper-surface mid-chord is sensitive to small shape changes (camber/
    // thickness pulse), so the derivative is large enough that central-diff
    // noise stays low relative to it.
    const std::size_t k = 0;
    const Eigen::MatrixXd& base_coords = base_stations[k].coords;
    const int N = static_cast<int>(base_coords.cols());
    const int j_node = N / 4;  // somewhere on the upper surface
    const int i_comp = 1;      // z component

    Scalar dCL_dcoord_analytic = g.dCL_dcoords[k](i_comp, j_node);
    Scalar dCD_dcoord_analytic = g.dCD_dcoords[k](i_comp, j_node);

    auto SolveWithPerturbedCoord = [&](Scalar h) -> CouplerResult {
        std::vector<CouplingStation> stations = BuildBaseStations();
        stations[k].coords(i_comp, j_node) += h;
        Aircraft a = BuildSmallWing(3);
        VortexLatticeSolver s;
        s.SetAircraft(a);
        FlightCondition f;
        f.alpha = alpha_base;
        ViscousVlmCoupler c(s, f, BuildOptions(stations));
        return c.Solve();
    };

    const Scalar h_coord = 1.0e-5;
    CouplerResult plus = SolveWithPerturbedCoord(+h_coord);
    CouplerResult minus = SolveWithPerturbedCoord(-h_coord);
    ASSERT_TRUE(plus.converged);
    ASSERT_TRUE(minus.converged);
    Scalar dCL_fd = (plus.results.CL - minus.results.CL) / (2.0 * h_coord);
    Scalar dCD_fd = (plus.results.CD - minus.results.CD) / (2.0 * h_coord);

    // Same 5e-6 floor as the α∞ test: dominated by the FD reference's
    // truncation/roundoff budget. dCL and dCD here are O(1) per coord
    // unit, so absolute tolerance ≈ 5e-6 too.
    Scalar cl_rel = std::abs(dCL_dcoord_analytic - dCL_fd) /
                    std::max(std::abs(dCL_fd), Scalar(1e-3));
    Scalar cd_rel = std::abs(dCD_dcoord_analytic - dCD_fd) /
                    std::max(std::abs(dCD_fd), Scalar(1e-4));
    // Shape-DV FD floor is dominated by Scirocco's internal primal-solve
    // precision (~1e-8 in Cl) divided by the coord FD step (1e-5):
    // central-diff noise ≈ 5e-4 relative. The analytical adjoint chain
    // itself is unlimited; this is the FD reference floor.
    EXPECT_LT(cl_rel, 5.0e-3)
        << "analytic dCL/dx = " << dCL_dcoord_analytic << ", FD = " << dCL_fd
        << ", rel err = " << cl_rel;
    EXPECT_LT(cd_rel, 5.0e-3)
        << "analytic dCD/dx = " << dCD_dcoord_analytic << ", FD = " << dCD_fd
        << ", rel err = " << cd_rel;
}

// Swept-wing primal: the local-frame α_eff projection should keep the
// coupler well-posed even with non-trivial sweep. NACA 0012 stations,
// uniform Re. The sweep angle here (tan⁻¹(2/5) ≈ 21.8°) is large enough
// to exercise the chord-plane projection.
TEST(ViscousCouplerTest, ConvergesOnSweptWing) {
    Aircraft aircraft = BuildSweptWing(/*n_spanwise=*/3, /*sweep_x=*/2.0);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition flight;
    flight.alpha = 2.0 * kDegreesToRadians;

    CouplerOptions opt;
    opt.stations = UniformStations({0.05, 0.5, 0.95}, /*panels=*/99);
    opt.mach = 0.0;
    opt.max_iterations = 20;
    opt.tolerance = 1.0e-6;

    ViscousVlmCoupler coupler(solver, flight, opt);
    CouplerResult res = coupler.Solve();

    EXPECT_TRUE(res.converged) << "Final residual " << res.final_residual
                               << " after " << res.iterations << " iterations";
    EXPECT_GT(res.results.CDv, 0.0);
    EXPECT_LT(res.results.CDv, 0.05);
    EXPECT_GT(res.results.CD, res.results.CDi_nearfield);
    EXPECT_GT(res.results.CL, 0.03)
        << "swept wing must produce real lift, CL = " << res.results.CL;
}

}  // namespace
}  // namespace pegasus::viscous
