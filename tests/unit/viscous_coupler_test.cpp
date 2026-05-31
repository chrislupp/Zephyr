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

namespace pegasus::viscous {
namespace {

// Build a small rectangular AR=10 wing with only enough strips to keep
// Scirocco-driven Newton iterations under a reasonable wall-clock budget.
// Each strip means another viscous airfoil solve per residual evaluation,
// and the FD-Jacobian builds N+1 evaluations per Newton step.
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

// Convenience: build a coupler around a freshly-meshed solver and return
// both. The coupler holds a reference to the solver, so the caller must
// keep `solver` alive at least as long as the returned coupler.
struct CouplerFixture {
    Aircraft aircraft;
    VortexLatticeSolver solver;
};

// Smoke test: the Newton loop converges on a small wing at low α and
// produces a viscous drag contribution that is positive and non-trivial.
TEST(ViscousCouplerTest, ConvergesOnSmallWingAndAddsViscousDrag) {
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/3);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition flight;
    flight.alpha = 2.0 * kDegreesToRadians;

    CouplerOptions opt;
    opt.airfoil_naca = "0012";
    opt.reynolds = 1.0e6;
    opt.mach = 0.0;
    opt.scirocco_panels = 99;  // smaller for faster test
    opt.max_iterations = 4;
    opt.tolerance = 5.0e-3;  // loose tolerance — smoke test, not a polar match

    ViscousVlmCoupler coupler(solver, flight, opt);
    CouplerResult res = coupler.Solve();

    EXPECT_TRUE(res.converged) << "Final residual " << res.final_residual
                               << " after " << res.iterations << " iterations";
    EXPECT_GT(res.results.CDv, 0.0);
    EXPECT_LT(res.results.CDv, 0.05);  // sanity bound

    // After viscous correction, total CD should exceed inviscid induced drag.
    EXPECT_GT(res.results.CD, res.results.CDi_nearfield);
}

// Validate the coupled total derivative dCD/dα∞ against a fully-converged
// finite difference over α∞. This exercises the implicit-function chain
// inside ComputeAlphaDerivatives(): J du/dα + ∂R/∂α = 0,  then
// dF/dα = ∂F/∂α + ∂F/∂u · du/dα.
TEST(ViscousCouplerTest, CoupledGradientMatchesFullFiniteDifference) {
    auto BuildOptions = []() {
        CouplerOptions opt;
        opt.airfoil_naca = "0012";
        opt.reynolds = 1.0e6;
        opt.mach = 0.0;
        opt.scirocco_panels = 99;
        opt.max_iterations = 5;
        opt.tolerance = 1.0e-4;
        // The inner FD step for ∂F/∂u and ∂R/∂α should be small enough to
        // give a clean directional derivative but large enough to swamp
        // Scirocco's own convergence noise; ~0.05° is a sweet spot.
        opt.fd_step_rad = 0.05 * kDegreesToRadians;
        return opt;
    };

    const Scalar alpha_base = 2.0 * kDegreesToRadians;
    const Scalar h_alpha = 0.25 * kDegreesToRadians;  // outer FD step
    CouplerOptions opt = BuildOptions();

    // --- Analytical (coupler implicit-function) gradient at α_base
    Aircraft aircraft = BuildSmallWing(/*n_spanwise=*/2);
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    FlightCondition flight;
    flight.alpha = alpha_base;

    ViscousVlmCoupler coupler(solver, flight, opt);
    CouplerResult base = coupler.Solve();
    ASSERT_TRUE(base.converged);

    CoupledGradients g = coupler.ComputeAlphaDerivatives();

    // --- Reference: re-converge the coupler at α_base ± h, take central
    // difference of CL and CD. Each independent re-solve gets its own
    // Aircraft/Solver to avoid carrying over the previously-mutated mesh
    // state.
    auto SolveAt = [&](Scalar alpha) -> CouplerResult {
        Aircraft a = BuildSmallWing(/*n_spanwise=*/2);
        VortexLatticeSolver s;
        s.SetAircraft(a);
        FlightCondition f;
        f.alpha = alpha;
        ViscousVlmCoupler c(s, f, BuildOptions());
        CouplerResult res = c.Solve();
        return res;
    };

    CouplerResult plus = SolveAt(alpha_base + h_alpha);
    CouplerResult minus = SolveAt(alpha_base - h_alpha);
    ASSERT_TRUE(plus.converged);
    ASSERT_TRUE(minus.converged);

    Scalar dCL_dalpha_fd =
        (plus.results.CL - minus.results.CL) / (2.0 * h_alpha);
    Scalar dCD_dalpha_fd =
        (plus.results.CD - minus.results.CD) / (2.0 * h_alpha);

    // 5% tolerance: FD-of-coupler vs implicit-function gradient agreement
    // is limited by the inner-FD Jacobian noise floor and the outer step's
    // truncation error.
    EXPECT_NEAR(g.dCL_dalpha, dCL_dalpha_fd,
                0.05 * std::max(std::abs(dCL_dalpha_fd), Scalar(1e-3)))
        << "analytic dCL/dα = " << g.dCL_dalpha
        << ", FD reference = " << dCL_dalpha_fd;
    EXPECT_NEAR(g.dCD_dalpha, dCD_dalpha_fd,
                0.10 * std::max(std::abs(dCD_dalpha_fd), Scalar(1e-4)))
        << "analytic dCD/dα = " << g.dCD_dalpha
        << ", FD reference = " << dCD_dalpha_fd;
}

}  // namespace
}  // namespace pegasus::viscous
