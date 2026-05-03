/*
    Project Pegasus
    Version 0.1.0

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

/// @file stability_analysis.cpp
/// @brief Example: Stability derivatives and dynamic mode analysis.

#include "pegasus/pegasus.h"

#include <iomanip>
#include <iostream>

int main() {
    using namespace pegasus;

    std::cout << "VLM Stability Analysis Example\n";
    std::cout << "==============================\n\n";

    // Create aircraft
    const double wing_span = 10.0;
    const double wing_chord = 1.5;
    const double wing_area = wing_span * wing_chord;

    // Wing with dihedral and sweep
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(16);
    wing.SetYDuplicate(0.0);

    Section wing_root(Vec3(0.0, 0.0, 0.0), wing_chord, 2.0);  // 2 deg twist
    Section wing_tip(Vec3(0.3, wing_span / 2.0, 0.3), wing_chord * 0.6, -1.0);
    wing.AddSection(wing_root);
    wing.AddSection(wing_tip);

    // Horizontal tail
    Surface htail("Horizontal Tail");
    htail.SetChordwisePanelCount(4);
    htail.SetSpanwisePanelCount(8);
    htail.SetYDuplicate(0.0);

    Section htail_root(Vec3(4.5, 0.0, 0.2), 0.6);
    Section htail_tip(Vec3(4.6, 1.5, 0.2), 0.4);
    htail.AddSection(htail_root);
    htail.AddSection(htail_tip);

    // Vertical tail
    Surface vtail("Vertical Tail");
    vtail.SetChordwisePanelCount(4);
    vtail.SetSpanwisePanelCount(6);

    Section vtail_root(Vec3(4.3, 0.0, 0.2), 0.8);
    Section vtail_tip(Vec3(4.6, 0.0, 1.2), 0.4);
    vtail.AddSection(vtail_root);
    vtail.AddSection(vtail_tip);

    Aircraft aircraft("Light Aircraft");
    aircraft.AddSurface(wing);
    aircraft.AddSurface(htail);
    aircraft.AddSurface(vtail);
    aircraft.SetReferenceArea(wing_area);
    aircraft.SetReferenceChord(wing_chord);
    aircraft.SetReferenceSpan(wing_span);
    aircraft.SetMomentReferencePoint(Vec3(0.4, 0.0, 0.0));

    std::cout << "Aircraft: " << aircraft.GetTitle() << "\n";
    std::cout << "  Reference area: " << wing_area << " m^2\n";
    std::cout << "  Reference chord: " << wing_chord << " m\n";
    std::cout << "  Reference span: " << wing_span << " m\n\n";

    // Solve and compute stability derivatives
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;
    condition.mach = 0.2;

    SolutionState solution = solver.Solve(condition);

    ForceCalculator force_calc;
    AerodynamicResults aero = force_calc.ComputeForces(
        aircraft, solver.GetMesh(), solution, condition);

    std::cout << "Flight Condition:\n";
    std::cout << "  Alpha: " << condition.alpha * kRadiansToDegrees << " deg\n";
    std::cout << "  Mach: " << condition.mach << "\n\n";

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Aerodynamic Coefficients:\n";
    std::cout << "  CL: " << aero.CL << "\n";
    std::cout << "  CD: " << aero.CD << "\n";
    std::cout << "  Cm: " << aero.Cm << "\n\n";

    // Compute stability derivatives
    StabilityDerivativeCalculator deriv_calc;
    StabilityDerivatives derivs =
        deriv_calc.ComputeDerivatives(aircraft, solver, condition);

    std::cout << std::setprecision(4);
    std::cout << "Longitudinal Stability Derivatives:\n";
    std::cout << "  CL_alpha: " << derivs.CL_alpha << " /rad\n";
    std::cout << "  Cm_alpha: " << derivs.Cm_alpha << " /rad\n";
    std::cout << "  Cm_q:     " << derivs.Cm_q << "\n\n";

    std::cout << "Lateral-Directional Stability Derivatives:\n";
    std::cout << "  CY_beta:  " << derivs.CY_beta << " /rad\n";
    std::cout << "  Cl_beta:  " << derivs.Cl_beta << " /rad (dihedral)\n";
    std::cout << "  Cn_beta:  " << derivs.Cn_beta << " /rad (directional)\n";
    std::cout << "  Cl_p:     " << derivs.Cl_p << " (roll damping)\n";
    std::cout << "  Cn_r:     " << derivs.Cn_r << " (yaw damping)\n\n";

    // Set up mass properties for eigenmode analysis
    MassProperties mass;
    mass.mass = 1000.0;
    mass.Ixx = 1000.0;
    mass.Iyy = 2000.0;
    mass.Izz = 2500.0;
    mass.Ixz = 100.0;

    FlightState flight;
    flight.velocity = 50.0;
    flight.alpha = condition.alpha;
    flight.density = 1.225;
    flight.gravity = 9.81;
    flight.theta = flight.alpha;

    std::cout << "Flight State for Dynamic Analysis:\n";
    std::cout << "  Velocity: " << flight.velocity << " m/s\n";
    std::cout << "  Dynamic pressure: " << std::setprecision(1)
              << flight.GetDynamicPressure() << " Pa\n\n";

    // Eigenmode analysis
    EigenmodeAnalysis eigenmode_analysis;

    EigenmodeResults lon_results =
        eigenmode_analysis.AnalyzeLongitudinal(aircraft, mass, derivs, flight);

    std::cout << "Longitudinal Dynamic Modes:\n";
    std::cout << "--------------------------\n";
    if (lon_results.is_valid) {
        for (const auto& mode : lon_results.modes) {
            std::cout << "  " << GetModeTypeName(mode.type) << ":\n";
            std::cout << std::setprecision(4);
            std::cout << "    Eigenvalue: " << mode.eigenvalue.real();
            if (mode.IsOscillatory()) {
                std::cout << " +/- " << std::abs(mode.eigenvalue.imag())
                          << "i\n";
                std::cout << "    Frequency: " << mode.GetFrequency()
                          << " rad/s\n";
                std::cout << "    Damping ratio: " << mode.GetDampingRatio()
                          << "\n";
            } else {
                std::cout << "\n";
            }
            std::cout << "    Stable: " << (mode.IsStable() ? "YES" : "NO")
                      << "\n\n";
        }
    }

    std::cout << "Stability Summary:\n";
    std::cout << "-----------------\n";
    std::cout << "  Static longitudinal: "
              << (derivs.Cm_alpha < 0 ? "STABLE" : "UNSTABLE") << "\n";
    std::cout << "  Directional: "
              << (derivs.Cn_beta > 0 ? "STABLE" : "UNSTABLE") << "\n";
    std::cout << "  Dihedral effect: "
              << (derivs.Cl_beta < 0 ? "Present" : "Absent") << "\n";

    std::cout << "\nStability analysis complete.\n";
    return 0;
}
