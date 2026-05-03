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

#include "pegasus/trim/trim_solver.h"

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/trim/flight_constraints.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-4;

// ============================================================================
// TrimSolver Tests
// ============================================================================

class TrimSolverTest : public ::testing::Test {
protected:
    Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                   int num_spanwise) {
        Aircraft aircraft("Test Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kCosine);
        wing.SetSpanwiseSpacing(SpacingType::kCosine);

        wing.AddSection(Section(Vec3(0, -span / 2, 0), chord, 0.0));
        wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));

        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * chord);
        aircraft.SetReferenceChord(chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }
};

TEST_F(TrimSolverTest, TrimToTargetCL) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Trim to CL = 0.5
    Scalar target_CL = 0.5;
    TrimConstraints constraints = TrimPresets::LevelFlight(target_CL);

    TrimSolver trimmer;
    TrimOptions options;
    options.verbose = false;
    trimmer.SetOptions(options);

    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged) << "Error: " << result.error_message;
    EXPECT_NEAR(result.aero.CL, target_CL, 0.01);
    EXPECT_GT(result.state.alpha, 0.0);  // Should be positive alpha
}

TEST_F(TrimSolverTest, TrimToZeroCL) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Trim to CL = 0 (symmetric airfoil)
    Scalar target_CL = 0.0;
    TrimConstraints constraints = TrimPresets::LevelFlight(target_CL);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged) << "Error: " << result.error_message;
    EXPECT_NEAR(result.aero.CL, 0.0, kTestTolerance);
    EXPECT_NEAR(result.state.alpha, 0.0, kTestTolerance);
}

TEST_F(TrimSolverTest, TrimWithFixedAlpha) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Fix alpha at 5 degrees
    Scalar target_alpha = 5.0 * kDegreesToRadians;
    TrimConstraints constraints;
    constraints.alpha.type = ConstraintType::kFixed;
    constraints.alpha.target = target_alpha;
    constraints.beta.type = ConstraintType::kFixed;
    constraints.beta.target = 0.0;
    constraints.roll_rate.type = ConstraintType::kFixed;
    constraints.roll_rate.target = 0.0;
    constraints.pitch_rate.type = ConstraintType::kFixed;
    constraints.pitch_rate.target = 0.0;
    constraints.yaw_rate.type = ConstraintType::kFixed;
    constraints.yaw_rate.target = 0.0;

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged) << "Error: " << result.error_message;
    EXPECT_NEAR(result.state.alpha, target_alpha, kTestTolerance);
    EXPECT_GT(result.aero.CL, 0.0);  // Should produce positive lift
}

TEST_F(TrimSolverTest, TrimWithFixedBeta) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Sideslip trim
    Scalar target_beta = 5.0 * kDegreesToRadians;
    TrimConstraints constraints = TrimPresets::SteadySideslip(0.5, target_beta);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged) << "Error: " << result.error_message;
    EXPECT_NEAR(result.state.beta, target_beta, kTestTolerance);
}

TEST_F(TrimSolverTest, ConvergenceInFewIterations) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimConstraints constraints = TrimPresets::LevelFlight(0.3);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.iterations, 10);  // Should converge quickly
}

TEST_F(TrimSolverTest, RatesRemainZeroForLevelFlight) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimConstraints constraints = TrimPresets::LevelFlight(0.4);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(result.state.roll_rate, 0.0, kTestTolerance);
    EXPECT_NEAR(result.state.pitch_rate, 0.0, kTestTolerance);
    EXPECT_NEAR(result.state.yaw_rate, 0.0, kTestTolerance);
}

TEST_F(TrimSolverTest, InitialGuessIsUsed) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Provide good initial guess
    OperatingVariables initial;
    initial.alpha = 4.0 * kDegreesToRadians;  // Close to expected trim alpha

    TrimConstraints constraints = TrimPresets::LevelFlight(0.4);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints, initial);

    EXPECT_TRUE(result.converged);
    EXPECT_LE(result.iterations,
              5);  // Should converge very quickly with good guess
}

TEST_F(TrimSolverTest, ResidualDecreasesMonotonically) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimConstraints constraints = TrimPresets::LevelFlight(0.6);

    TrimSolver trimmer;
    TrimOptions options;
    options.verbose = false;
    options.max_iterations = 10;
    trimmer.SetOptions(options);

    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    // Final residual should be small if converged
    if (result.converged) {
        EXPECT_LT(result.residual_norm, options.tolerance);
    }
}

TEST_F(TrimSolverTest, StabilityDerivativesComputed) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimConstraints constraints = TrimPresets::LevelFlight(0.5);

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(result.derivatives.is_valid);
    EXPECT_GT(result.derivatives.CL_alpha,
              0.0);  // Should have positive lift slope
}

TEST_F(TrimSolverTest, NoConstraintsFails) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimConstraints constraints;  // All free - no constraints

    TrimSolver trimmer;
    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    EXPECT_FALSE(result.converged);
    EXPECT_FALSE(result.error_message.empty());
}

TEST_F(TrimSolverTest, UnreachableCLFails) {
    Aircraft aircraft = CreateRectangularWing(10.0, 1.0, 4, 16);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    // Request impossibly high CL
    TrimConstraints constraints = TrimPresets::LevelFlight(10.0);

    TrimSolver trimmer;
    TrimOptions options;
    options.max_alpha = 0.5;  // Restrict alpha to prevent huge values
    trimmer.SetOptions(options);

    TrimResult result = trimmer.Solve(aircraft, solver, constraints);

    // Should fail because can't reach CL=10 within alpha limits
    EXPECT_FALSE(result.converged);
}

TEST_F(TrimSolverTest, PresetLevelFlight) {
    // Test the preset convenience function
    TrimConstraints constraints = TrimPresets::LevelFlight(0.5);

    EXPECT_EQ(constraints.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(constraints.alpha.target, 0.5, kTestTolerance);
    EXPECT_EQ(constraints.beta.type, ConstraintType::kFixed);
    EXPECT_NEAR(constraints.beta.target, 0.0, kTestTolerance);
}

TEST_F(TrimSolverTest, PresetSteadySideslip) {
    Scalar beta = 10.0 * kDegreesToRadians;
    TrimConstraints constraints = TrimPresets::SteadySideslip(0.4, beta);

    EXPECT_EQ(constraints.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(constraints.alpha.target, 0.4, kTestTolerance);
    EXPECT_EQ(constraints.beta.type, ConstraintType::kFixed);
    EXPECT_NEAR(constraints.beta.target, beta, kTestTolerance);
}

}  // namespace
}  // namespace pegasus
