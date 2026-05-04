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

#include "pegasus/geometry/mesh_input.h"

#include "pegasus/aero/force_calculator.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

/// Build a flat rectangular wing mesh with given dimensions.
/// LE at x=0, TE at x=chord. Span from y=0 to y=span.
SurfaceMesh BuildRectangularMesh(Index nx, Index ny, Scalar chord,
                                 Scalar span) {
    SurfaceMesh sm;
    sm.name = "Wing";
    sm.nx = nx;
    sm.ny = ny;
    sm.nodes.resize((nx + 1) * (ny + 1), 3);

    for (Index i = 0; i <= nx; ++i) {
        for (Index j = 0; j <= ny; ++j) {
            Scalar x = chord * static_cast<Scalar>(i) / static_cast<Scalar>(nx);
            Scalar y = span * static_cast<Scalar>(j) / static_cast<Scalar>(ny);
            Scalar z = 0.0;
            sm.SetNode(i, j, Vec3(x, y, z));
        }
    }
    return sm;
}

/// Build an Aircraft object with a rectangular wing matching the mesh.
Aircraft BuildRectangularAircraft(Index nx, Index ny, Scalar chord,
                                  Scalar span) {
    Aircraft aircraft("RectWing");
    aircraft.SetReferenceArea(chord * span);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);
    aircraft.SetMomentReferencePoint(Vec3(0.25 * chord, 0.5 * span, 0));

    Surface surface("Wing");
    surface.SetChordwisePanelCount(nx);
    surface.SetChordwiseSpacing(SpacingType::kUniform);
    surface.SetSpanwisePanelCount(ny);
    surface.SetSpanwiseSpacing(SpacingType::kUniform);

    Section sec_root;
    sec_root.SetLeadingEdgePosition(Vec3(0, 0, 0));
    sec_root.SetChord(chord);
    sec_root.SetTwistAngle(0.0);
    surface.AddSection(sec_root);

    Section sec_tip;
    sec_tip.SetLeadingEdgePosition(Vec3(0, span, 0));
    sec_tip.SetChord(chord);
    sec_tip.SetTwistAngle(0.0);
    surface.AddSection(sec_tip);

    aircraft.AddSurface(surface);
    return aircraft;
}

// Test that mesh-based solver produces valid results
TEST(MeshInputTest, BasicSolveProducesValidResults) {
    const Index nx = 4;
    const Index ny = 8;
    const Scalar chord = 1.0;
    const Scalar span = 4.0;

    SurfaceMesh sm = BuildRectangularMesh(nx, ny, chord, span);

    MeshInput input;
    input.surfaces.push_back(sm);
    input.reference.reference_area = chord * span;
    input.reference.reference_chord = chord;
    input.reference.reference_span = span;
    input.reference.moment_reference_point = Vec3(0.25 * chord, 0.5 * span, 0);

    VortexLatticeSolver solver;
    solver.SetMesh(input);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState state = solver.Solve(condition);
    ASSERT_TRUE(state.is_valid);
    EXPECT_GT(state.circulation.size(), 0);

    ForceCalculator calc;
    AerodynamicResults results =
        calc.ComputeForces(input.reference, solver.GetMesh(), state, condition);

    EXPECT_TRUE(results.is_valid);
    // A flat plate at 5 degrees should have positive lift
    EXPECT_GT(results.CL, 0.0);
    // CL should be reasonable (thin-airfoil theory ~2*pi*alpha ~ 0.55)
    EXPECT_GT(results.CL, 0.2);
    EXPECT_LT(results.CL, 1.0);
}

// Test that mesh path and aircraft path give similar results
TEST(MeshInputTest, ConsistencyWithAircraftPath) {
    const Index nx = 6;
    const Index ny = 12;
    const Scalar chord = 1.0;
    const Scalar span = 6.0;

    // Mesh-based path
    SurfaceMesh sm = BuildRectangularMesh(nx, ny, chord, span);
    MeshInput input;
    input.surfaces.push_back(sm);
    input.reference.reference_area = chord * span;
    input.reference.reference_chord = chord;
    input.reference.reference_span = span;
    input.reference.moment_reference_point = Vec3(0.25, 3.0, 0);

    VortexLatticeSolver solver_mesh;
    solver_mesh.SetMesh(input);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState state_mesh = solver_mesh.Solve(condition);
    ForceCalculator calc;
    AerodynamicResults results_mesh = calc.ComputeForces(
        input.reference, solver_mesh.GetMesh(), state_mesh, condition);

    // Aircraft-based path
    Aircraft aircraft = BuildRectangularAircraft(nx, ny, chord, span);
    VortexLatticeSolver solver_aircraft;
    solver_aircraft.SetAircraft(aircraft);

    SolutionState state_aircraft = solver_aircraft.Solve(condition);
    AerodynamicResults results_aircraft = calc.ComputeForces(
        aircraft, solver_aircraft.GetMesh(), state_aircraft, condition);

    // Both should give positive, valid results
    ASSERT_TRUE(results_mesh.is_valid);
    ASSERT_TRUE(results_aircraft.is_valid);

    // CL should be reasonably close. Note: the mesh-based path uses uniform
    // spacing directly, while the aircraft path may differ slightly in how
    // it handles the 1/4-3/4 chord placement. But for uniform spacing
    // and a flat plate, they should be very close.
    EXPECT_NEAR(results_mesh.CL, results_aircraft.CL, 0.05);
    EXPECT_NEAR(results_mesh.Cm, results_aircraft.Cm, 0.05);
}

// Test mesh generation dimensions
TEST(MeshInputTest, MeshDimensions) {
    const Index nx = 3;
    const Index ny = 5;

    SurfaceMesh sm = BuildRectangularMesh(nx, ny, 1.0, 4.0);

    MeshInput input;
    input.surfaces.push_back(sm);

    VortexLatticeSolver solver;
    solver.SetMesh(input);

    const VortexLattice& mesh = solver.GetMesh();
    EXPECT_EQ(mesh.GetNumVortices(), nx * ny);
    EXPECT_EQ(mesh.GetNumStrips(), ny);
    EXPECT_EQ(mesh.GetNumSurfaces(), 1);
    EXPECT_EQ(mesh.GetSurfaceInfo(0).num_chordwise, nx);
    EXPECT_EQ(mesh.GetSurfaceInfo(0).num_spanwise, ny);
}

}  // namespace
}  // namespace pegasus
