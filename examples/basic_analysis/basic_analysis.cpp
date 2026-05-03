/*
    Project Pegasus
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

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

/// @file basic_analysis.cpp
/// @brief Example: Basic aerodynamic analysis of a rectangular wing.

#include "pegasus/pegasus.h"

#include <iomanip>
#include <iostream>

int main() {
    using namespace pegasus;

    std::cout << "VLM Basic Analysis Example\n";
    std::cout << "==========================\n\n";

    // Wing parameters
    const double span = 10.0;
    const double chord = 1.0;
    const double aspect_ratio = span / chord;

    std::cout << "Creating rectangular wing:\n";
    std::cout << "  Span: " << span << " m\n";
    std::cout << "  Chord: " << chord << " m\n";
    std::cout << "  Aspect ratio: " << aspect_ratio << "\n\n";

    // Create wing surface
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(20);
    wing.SetChordwiseSpacing(SpacingType::kCosine);
    wing.SetSpanwiseSpacing(SpacingType::kCosine);
    wing.SetYDuplicate(0.0);

    Section root(Vec3(0.0, 0.0, 0.0), chord);
    Section tip(Vec3(0.0, span / 2.0, 0.0), chord);
    wing.AddSection(root);
    wing.AddSection(tip);

    // Create aircraft
    Aircraft aircraft("Rectangular Wing");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(span * chord);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);

    // Set up solver
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    std::cout << "Vortex lattice generated:\n";
    std::cout << "  Vortex elements: " << solver.GetMesh().GetNumVortices()
              << "\n";
    std::cout << "  Strips: " << solver.GetMesh().GetNumStrips() << "\n\n";

    // Solve at multiple angles of attack
    std::cout << "Angle of Attack Sweep\n";
    std::cout << "---------------------\n";
    std::cout << std::setw(10) << "Alpha" << std::setw(12) << "CL"
              << std::setw(12) << "CD" << std::setw(12) << "L/D" << "\n";

    ForceCalculator calculator;

    for (double alpha_deg = -2.0; alpha_deg <= 10.0; alpha_deg += 2.0) {
        FlightCondition condition;
        condition.alpha = alpha_deg * kDegreesToRadians;

        SolutionState solution = solver.Solve(condition);
        AerodynamicResults results = calculator.ComputeForces(
            aircraft, solver.GetMesh(), solution, condition);

        double ld_ratio = (results.CD > 1e-6) ? results.CL / results.CD : 0.0;

        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::setw(10) << alpha_deg;
        std::cout << std::setprecision(5);
        std::cout << std::setw(12) << results.CL << std::setw(12) << results.CD
                  << std::setw(12) << ld_ratio << "\n";
    }

    // Comparison with theory
    std::cout << "\nComparison with Theory\n";
    std::cout << "----------------------\n";

    FlightCondition cond;
    cond.alpha = 5.0 * kDegreesToRadians;
    SolutionState sol = solver.Solve(cond);
    AerodynamicResults res =
        calculator.ComputeForces(aircraft, solver.GetMesh(), sol, cond);

    double CL_alpha_2D = 2.0 * kPi;
    double e = 0.95;
    double CL_alpha_theory =
        CL_alpha_2D / (1.0 + CL_alpha_2D / (kPi * e * aspect_ratio));
    double CL_alpha_computed = res.CL / cond.alpha;

    std::cout << std::setprecision(4);
    std::cout << "  Theoretical CL_alpha: " << CL_alpha_theory << " /rad\n";
    std::cout << "  Computed CL_alpha:    " << CL_alpha_computed << " /rad\n";
    std::cout << "  Difference: " << std::setprecision(2)
              << 100.0 * (CL_alpha_computed - CL_alpha_theory) / CL_alpha_theory
              << "%\n";

    std::cout << "\nAnalysis complete.\n";
    return 0;
}
