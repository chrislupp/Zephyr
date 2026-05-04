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

/// @file vtk_visualization.cpp
/// @brief Example: Export VLM results to VTK format for visualization in
/// ParaView

#include "pegasus/aero/force_calculator.h"
#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/visualization/vtk_export.h"

#include <filesystem>
#include <iostream>

int main(int argc, char* argv[]) {
    using namespace pegasus;

    std::cout << "VLM VTK Visualization Export Example\n";
    std::cout << "=====================================\n\n";

    // Check VTK availability
    if (!VtkExporter::IsAvailable()) {
        std::cerr << "Error: VTK support not available.\n";
        std::cerr << "Rebuild with VTK: cmake -DVLM_ENABLE_VTK=ON ..\n";
        std::cerr << "Make sure VTK is installed on your system.\n";
        return 1;
    }

    std::cout << "VTK support: Available\n\n";

    // Output directory
    std::string output_dir = (argc > 1) ? argv[1] : "./vtk_output";
    std::cout << "Output directory: " << output_dir << "\n\n";

    // Create a tapered wing with washout
    std::cout << "Creating aircraft geometry...\n";

    Surface wing("Wing");
    wing.SetChordwisePanelCount(8);
    wing.SetSpanwisePanelCount(16);
    wing.SetYDuplicate(0.0);  // Mirror about Y=0

    // Root section: chord=1.0, incidence=2 deg
    Section root(Vec3(0.0, 0.0, 0.0), 1.0, 2.0 * kDegreesToRadians);
    // Tip section: chord=0.5, incidence=-1 deg (washout)
    Section tip(Vec3(0.25, 5.0, 0.0), 0.5, -1.0 * kDegreesToRadians);

    wing.AddSection(root);
    wing.AddSection(tip);

    Aircraft aircraft("Tapered Wing with Washout");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(7.5);    // Approximate planform area
    aircraft.SetReferenceChord(0.75);  // Mean chord
    aircraft.SetReferenceSpan(10.0);   // Full span

    std::cout << "  Reference area: " << aircraft.GetReferenceArea()
              << " m^2\n";
    std::cout << "  Reference chord: " << aircraft.GetReferenceChord()
              << " m\n";
    std::cout << "  Reference span: " << aircraft.GetReferenceSpan()
              << " m\n\n";

    // Set up solver
    std::cout << "Setting up VLM solver...\n";
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    std::cout << "  Mesh vortices: " << solver.GetMesh().GetNumVortices()
              << "\n";
    std::cout << "  Mesh strips: " << solver.GetMesh().GetNumStrips() << "\n\n";

    // Flight condition
    FlightCondition condition;
    condition.alpha = 6.0 * kDegreesToRadians;
    condition.beta = 0.0;

    std::cout << "Solving at:\n";
    std::cout << "  Alpha: " << condition.alpha * kRadiansToDegrees << " deg\n";
    std::cout << "  Beta: " << condition.beta * kRadiansToDegrees << " deg\n\n";

    // Solve
    SolutionState solution = solver.Solve(condition);
    if (!solution.is_valid) {
        std::cerr << "Error: VLM solution failed!\n";
        return 1;
    }

    // Compute forces
    ForceCalculator calculator;
    AerodynamicResults results = calculator.ComputeForces(
        aircraft, solver.GetMesh(), solution, condition);

    std::cout << "Aerodynamic Results:\n";
    std::cout << "  CL = " << results.CL << "\n";
    std::cout << "  CD = " << results.CD << " (near-field)\n";
    std::cout << "  CDi = " << results.CDi_ff << " (far-field induced)\n";
    std::cout << "  L/D = " << results.CL / results.CDi_ff << "\n";
    std::cout << "  Cm = " << results.Cm << "\n";
    std::cout << "  Span efficiency: " << results.span_efficiency << "\n\n";

    // Export to VTK
    std::cout << "Exporting VTK files...\n";

    VtkExportOptions options;
    options.export_cp = true;
    options.export_circulation = true;
    options.export_normals = true;
    options.export_indices = true;
    options.export_wake = true;
    options.wake_length_chords = 8.0;
    options.export_induced_velocity = true;
    options.export_spanwise_distributions = true;
    options.binary_format = true;  // Compact binary format

    bool success = VtkExporter::ExportAll(
        output_dir, solver.GetMesh(), solution, results, condition,
        aircraft.GetReferenceChord(), options);

    if (success) {
        std::cout << "\nExport successful!\n\n";
        std::cout << "Files created:\n";
        std::cout << "  - " << output_dir << "/surface_mesh.vtu\n";
        std::cout << "      Panel mesh colored by pressure coefficient (Cp)\n";
        std::cout
            << "      Cell data: Cp, Circulation, Surface, Strip, Normal\n\n";
        std::cout << "  - " << output_dir << "/spanwise_loads.vtp\n";
        std::cout << "      Spanwise distributions (CL, CD, moments)\n";
        std::cout << "      Point data: CL, CD, CY, SpanLoading, InducedAlpha, "
                     "CM, InducedVelocity\n\n";
        std::cout << "  - " << output_dir << "/wake.vtp\n";
        std::cout << "      Trailing vortex wake filaments\n";
        std::cout << "      Cell data: Circulation, Surface\n\n";
        std::cout << "  - " << output_dir << "/velocity_field.vtp\n";
        std::cout << "      Induced velocity vectors at control points\n";
        std::cout << "      Point data: InducedVelocity, TotalVelocity, "
                     "VelocityMagnitude\n\n";
        std::cout << "To visualize in ParaView:\n";
        std::cout << "  paraview " << output_dir << "/*.vt*\n\n";
        std::cout << "ParaView tips:\n";
        std::cout
            << "  - Use 'Glyph' filter on velocity_field.vtp to show arrows\n";
        std::cout
            << "  - Color surface_mesh.vtu by 'Cp' for pressure distribution\n";
        std::cout << "  - Use 'Tube' filter on wake.vtp for better wake "
                     "visualization\n";
    } else {
        std::cerr << "\nExport failed!\n";
        return 1;
    }

    return 0;
}
