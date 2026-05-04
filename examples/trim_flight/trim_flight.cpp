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

/// @file trim_flight.cpp
/// @brief Example: Trim solver for level flight conditions.

#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/trim/flight_constraints.h"
#include "pegasus/trim/trim_solver.h"

#include <iomanip>
#include <iostream>

int main() {
    using namespace pegasus;

    std::cout << "VLM Trim Flight Example\n";
    std::cout << "=======================\n\n";

    // Create aircraft with wing and tail
    const double wing_span = 10.0;
    const double wing_chord = 1.5;
    const double wing_area = wing_span * wing_chord;
    const double tail_x = 4.0;

    // Main wing
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(16);
    wing.SetYDuplicate(0.0);

    Section wing_root(Vec3(0.0, 0.0, 0.0), wing_chord);
    Section wing_tip(Vec3(0.0, wing_span / 2.0, 0.0), wing_chord);
    wing.AddSection(wing_root);
    wing.AddSection(wing_tip);

    // Horizontal tail
    Surface htail("Horizontal Tail");
    htail.SetChordwisePanelCount(4);
    htail.SetSpanwisePanelCount(8);
    htail.SetYDuplicate(0.0);

    Section tail_root(Vec3(tail_x, 0.0, 0.0), 0.6);
    Section tail_tip(Vec3(tail_x, 1.5, 0.0), 0.4);
    htail.AddSection(tail_root);
    htail.AddSection(tail_tip);

    Aircraft aircraft("Simple Aircraft");
    aircraft.AddSurface(wing);
    aircraft.AddSurface(htail);
    aircraft.SetReferenceArea(wing_area);
    aircraft.SetReferenceChord(wing_chord);
    aircraft.SetReferenceSpan(wing_span);
    aircraft.SetMomentReferencePoint(Vec3(0.25 * wing_chord, 0.0, 0.0));

    std::cout << "Aircraft Configuration:\n";
    std::cout << "  Wing span: " << wing_span << " m\n";
    std::cout << "  Wing chord: " << wing_chord << " m\n";
    std::cout << "  Tail arm: " << tail_x << " m\n\n";

    // Set up solver
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    std::cout << "Total vortex elements: " << solver.GetMesh().GetNumVortices()
              << "\n\n";

    // Trim to different CL targets
    TrimSolver trim_solver;

    std::vector<double> target_CLs = {0.3, 0.5, 0.7, 1.0};

    std::cout << "Trim Results for Level Flight\n";
    std::cout << "-----------------------------\n";
    std::cout << std::setw(10) << "CL_target" << std::setw(12) << "Alpha (deg)"
              << std::setw(12) << "CL" << std::setw(12) << "Cm" << std::setw(12)
              << "Iterations" << "\n";

    for (double target_CL : target_CLs) {
        TrimConstraints constraints = TrimPresets::LevelFlight(target_CL);
        OperatingVariables initial;
        initial.alpha = target_CL / 5.0;

        TrimResult result =
            trim_solver.Solve(aircraft, solver, constraints, initial);

        if (result.converged) {
            std::cout << std::fixed << std::setprecision(3);
            std::cout << std::setw(10) << target_CL;
            std::cout << std::setw(12)
                      << result.state.alpha * kRadiansToDegrees;
            std::cout << std::setprecision(5);
            std::cout << std::setw(12) << result.aero.CL;
            std::cout << std::setw(12) << result.aero.Cm;
            std::cout << std::setw(12) << result.iterations << "\n";
        } else {
            std::cout << std::setw(10) << target_CL << "  FAILED\n";
        }
    }

    std::cout << "\nTrim analysis complete.\n";
    return 0;
}
