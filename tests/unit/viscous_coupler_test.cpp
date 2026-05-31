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

}  // namespace
}  // namespace pegasus::viscous
