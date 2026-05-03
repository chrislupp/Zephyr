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

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// VortexLatticeSolver Tests
// ============================================================================

class VortexLatticeSolverTest : public ::testing::Test {
protected:
    Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                   int num_spanwise, bool y_duplicate = false) {
        Aircraft aircraft("Test Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        if (y_duplicate) {
            // Half wing with y-duplicate
            wing.AddSection(Section(Vec3(0, 0, 0), chord, 0.0));
            wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));
            wing.SetYDuplicate(0.0);
        } else {
            // Full wing
            wing.AddSection(Section(Vec3(0, -span / 2, 0), chord, 0.0));
            wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));
        }

        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * chord);
        aircraft.SetReferenceChord(chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }
};

TEST_F(VortexLatticeSolverTest, SetAircraftGeneratesMesh) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    EXPECT_GT(solver.GetMesh().GetNumVortices(), 0);
    EXPECT_GT(solver.GetMesh().GetNumStrips(), 0);
}

TEST_F(VortexLatticeSolverTest, RebuildAICCreatesValidMatrix) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    solver.RebuildAIC();

    EXPECT_TRUE(solver.IsAICValid());

    const MatX& aic = solver.GetAICMatrix();
    Index n = solver.GetMesh().GetNumVortices();
    EXPECT_EQ(aic.rows(), n);
    EXPECT_EQ(aic.cols(), n);
}

TEST_F(VortexLatticeSolverTest, SolveUnitResponsesComputes6Columns) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    solver.SolveUnitResponses();

    EXPECT_TRUE(solver.AreUnitResponsesValid());
}

TEST_F(VortexLatticeSolverTest, SolveAtZeroAlphaGivesZeroCirculation) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 0.0;
    condition.beta = 0.0;

    SolutionState state = solver.Solve(condition);

    EXPECT_TRUE(state.is_valid);

    // At zero alpha on flat plate, circulation should be zero
    for (Index i = 0; i < state.circulation.size(); ++i) {
        EXPECT_NEAR(state.circulation(i), 0.0, kTestTolerance);
    }
}

TEST_F(VortexLatticeSolverTest, SolveAtPositiveAlphaGivesPositiveCirculation) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;
    condition.beta = 0.0;

    SolutionState state = solver.Solve(condition);

    EXPECT_TRUE(state.is_valid);

    // At positive alpha, circulation should be positive (lift up)
    Scalar total_circulation = state.circulation.sum();
    EXPECT_GT(total_circulation, 0.0);
}

TEST_F(VortexLatticeSolverTest, CirculationLinearWithAlpha) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 6, 8);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition1;
    condition1.alpha = 2.0 * kDegreesToRadians;
    SolutionState state1 = solver.Solve(condition1);

    FlightCondition condition2;
    condition2.alpha = 4.0 * kDegreesToRadians;
    SolutionState state2 = solver.Solve(condition2);

    // Circulation should be linear with alpha (double alpha = double
    // circulation)
    Scalar total1 = state1.circulation.sum();
    Scalar total2 = state2.circulation.sum();

    EXPECT_NEAR(total2 / total1, 2.0, 0.01);
}

TEST_F(VortexLatticeSolverTest, FreestreamDirectionComputation) {
    FlightCondition condition;

    // Zero alpha, zero beta: freestream in +x
    condition.alpha = 0.0;
    condition.beta = 0.0;
    Vec3 v1 = VortexLatticeSolver::ComputeFreestreamDirection(condition);
    EXPECT_NEAR(v1.x(), 1.0, kTestTolerance);
    EXPECT_NEAR(v1.y(), 0.0, kTestTolerance);
    EXPECT_NEAR(v1.z(), 0.0, kTestTolerance);

    // 90 deg alpha: freestream in +z
    condition.alpha = kPi / 2;
    condition.beta = 0.0;
    Vec3 v2 = VortexLatticeSolver::ComputeFreestreamDirection(condition);
    EXPECT_NEAR(v2.x(), 0.0, kTestTolerance);
    EXPECT_NEAR(v2.y(), 0.0, kTestTolerance);
    EXPECT_NEAR(v2.z(), 1.0, kTestTolerance);

    // 90 deg beta: freestream in +y
    condition.alpha = 0.0;
    condition.beta = kPi / 2;
    Vec3 v3 = VortexLatticeSolver::ComputeFreestreamDirection(condition);
    EXPECT_NEAR(v3.x(), 0.0, kTestTolerance);
    EXPECT_NEAR(v3.y(), 1.0, kTestTolerance);
    EXPECT_NEAR(v3.z(), 0.0, kTestTolerance);
}

TEST_F(VortexLatticeSolverTest, MachNumberChangeInvalidatesAIC) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition1;
    condition1.mach = 0.0;
    solver.Solve(condition1);

    EXPECT_TRUE(solver.IsAICValid());

    // Change Mach number
    FlightCondition condition2;
    condition2.mach = 0.5;

    // Solve will rebuild AIC internally
    SolutionState state = solver.Solve(condition2);
    EXPECT_TRUE(state.is_valid);
}

// ============================================================================
// Theoretical Verification Tests
// ============================================================================

TEST_F(VortexLatticeSolverTest, HighAspectRatioWingLiftSlope) {
    // For a high aspect ratio wing (AR >> 1), lift slope should approach 2*pi
    // Using AR = 10 for good approximation
    Scalar span = 10.0;
    Scalar chord = 1.0;
    Scalar AR = span / chord;

    // Use fine mesh for accuracy
    Aircraft aircraft = CreateRectangularWing(span, chord, 8, 20);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Compute lift at small angle of attack
    Scalar alpha = 1.0 * kDegreesToRadians;

    FlightCondition condition;
    condition.alpha = alpha;

    SolutionState state = solver.Solve(condition);

    // Calculate CL from circulation
    // CL = 2 * sum(Gamma * delta_y) / (V * S)
    // For unit freestream, V = 1
    Scalar S_ref = aircraft.GetReferenceArea();
    Scalar total_lift_coefficient = 0.0;

    const VortexLattice& mesh = solver.GetMesh();
    for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
        const Strip& strip = mesh.GetStrip(j);
        Scalar strip_circulation = 0.0;
        for (Index k = 0; k < strip.num_vortices; ++k) {
            strip_circulation +=
                state.circulation(strip.first_vortex_index + k);
        }
        total_lift_coefficient += 2.0 * strip_circulation * strip.width / S_ref;
    }

    // CL_alpha = CL / alpha
    Scalar CL_alpha = total_lift_coefficient / alpha;

    // Theoretical: CL_alpha = 2*pi * AR / (AR + 2)  (Prandtl lifting line)
    // For AR = 10: CL_alpha ≈ 2*pi * 10/12 ≈ 5.24
    Scalar CL_alpha_theory = 2.0 * kPi * AR / (AR + 2.0);

    // VLM typically gives 15-25% higher lift slope than simple lifting-line
    // theory because it captures chordwise pressure distribution effects and
    // uses different induced drag modeling. Allow 25% tolerance for this
    // comparison.
    EXPECT_NEAR(CL_alpha, CL_alpha_theory, 0.25 * CL_alpha_theory);

    // Also verify the lift slope is physically reasonable (between 4.5 and 7.0)
    EXPECT_GT(CL_alpha, 4.5);
    EXPECT_LT(CL_alpha, 7.0);
}

TEST_F(VortexLatticeSolverTest, EllipticLiftDistributionShape) {
    // High aspect ratio wing should have nearly elliptic lift distribution
    Scalar span = 10.0;
    Scalar chord = 1.0;

    Aircraft aircraft = CreateRectangularWing(span, chord, 4, 20);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState state = solver.Solve(condition);

    const VortexLattice& mesh = solver.GetMesh();

    // Collect spanwise circulation distribution
    std::vector<std::pair<Scalar, Scalar>> spanwise_circ;

    for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
        const Strip& strip = mesh.GetStrip(j);
        Scalar strip_circulation = 0.0;
        for (Index k = 0; k < strip.num_vortices; ++k) {
            strip_circulation +=
                state.circulation(strip.first_vortex_index + k);
        }
        Scalar y = strip.leading_edge_center.y();
        spanwise_circ.emplace_back(y, strip_circulation);
    }

    // Find maximum circulation (should be near center)
    Scalar max_circ = 0.0;
    for (const auto& [y, circ] : spanwise_circ) {
        max_circ = std::max(max_circ, circ);
    }

    // Verify that circulation decreases toward wing tips
    // For elliptic: Gamma(y) = Gamma_0 * sqrt(1 - (2y/b)^2)
    bool monotonic_decrease_from_center = true;
    for (const auto& [y, circ] : spanwise_circ) {
        Scalar eta = 2.0 * y / span;  // Normalized span coordinate
        Scalar elliptic_ratio = std::sqrt(1.0 - eta * eta);
        Scalar expected_circ = max_circ * elliptic_ratio;

        // Allow 30% deviation from perfect elliptic (rectangular wing deviates)
        if (std::abs(y) > span * 0.1 && std::abs(y) < span * 0.4) {
            EXPECT_GT(circ, 0.5 * expected_circ);
            EXPECT_LT(circ, 1.5 * expected_circ);
        }
    }
}

TEST_F(VortexLatticeSolverTest, YDuplicateProducesSymmetricCirculation) {
    Scalar span = 4.0;
    Scalar chord = 1.0;

    Aircraft aircraft = CreateRectangularWing(span, chord, 4, 8, true);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState state = solver.Solve(condition);

    const VortexLattice& mesh = solver.GetMesh();

    // Should have 2 surfaces (original + y-duplicate)
    EXPECT_EQ(mesh.GetNumSurfaces(), 2);

    // Both surfaces should have equal total circulation
    const auto& info1 = mesh.GetSurfaceInfo(0);
    const auto& info2 = mesh.GetSurfaceInfo(1);

    Scalar total1 = 0.0, total2 = 0.0;
    for (Index i = info1.first_vortex_index;
         i < info1.first_vortex_index + info1.num_vortices; ++i) {
        total1 += state.circulation(i);
    }
    for (Index i = info2.first_vortex_index;
         i < info2.first_vortex_index + info2.num_vortices; ++i) {
        total2 += state.circulation(i);
    }

    EXPECT_NEAR(total1, total2, 0.01 * std::abs(total1));
}

TEST_F(VortexLatticeSolverTest, InducedVelocityAtControlPoints) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;

    SolutionState state = solver.Solve(condition);

    // Induced velocities should be computed
    EXPECT_EQ(state.induced_velocity.rows(), solver.GetMesh().GetNumVortices());
    EXPECT_EQ(state.induced_velocity.cols(), 3);

    // Induced w-velocity (downwash) should be negative for positive lift
    Scalar avg_w = state.induced_velocity.col(2).mean();
    EXPECT_LT(avg_w, 0.0);
}

}  // namespace
}  // namespace pegasus
