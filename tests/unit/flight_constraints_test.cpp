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

#include "pegasus/trim/flight_constraints.h"

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTol = 1.0e-10;

// ============================================================================
// OperatingVariables Tests
// ============================================================================

TEST(OperatingVariablesTest, DefaultNumVariables) {
    OperatingVariables vars;
    EXPECT_EQ(vars.GetNumVariables(), 5);
}

TEST(OperatingVariablesTest, NumVariablesWithControls) {
    OperatingVariables vars;
    vars.control_deflections.emplace_back("aileron", 0.0);
    EXPECT_EQ(vars.GetNumVariables(), 6);

    vars.control_deflections.emplace_back("elevator", 0.1);
    EXPECT_EQ(vars.GetNumVariables(), 7);
}

TEST(OperatingVariablesTest, DefaultValues) {
    OperatingVariables vars;
    EXPECT_NEAR(vars.alpha, 0.0, kTol);
    EXPECT_NEAR(vars.beta, 0.0, kTol);
    EXPECT_NEAR(vars.roll_rate, 0.0, kTol);
    EXPECT_NEAR(vars.pitch_rate, 0.0, kTol);
    EXPECT_NEAR(vars.yaw_rate, 0.0, kTol);
    EXPECT_TRUE(vars.control_deflections.empty());
}

// ============================================================================
// TrimConstraints Tests
// ============================================================================

TEST(TrimConstraintsTest, DefaultAllFree) {
    TrimConstraints constraints;
    EXPECT_EQ(constraints.alpha.type, ConstraintType::kFree);
    EXPECT_EQ(constraints.beta.type, ConstraintType::kFree);
    EXPECT_EQ(constraints.roll_rate.type, ConstraintType::kFree);
    EXPECT_EQ(constraints.pitch_rate.type, ConstraintType::kFree);
    EXPECT_EQ(constraints.yaw_rate.type, ConstraintType::kFree);
    EXPECT_EQ(constraints.GetNumConstraints(), 0);
}

TEST(TrimConstraintsTest, ConstraintCountWithControls) {
    TrimConstraints constraints;
    constraints.alpha.type = ConstraintType::kCL;
    constraints.alpha.target = 0.5;

    VariableConstraint ctrl;
    ctrl.type = ConstraintType::kPitchMoment;
    ctrl.target = 0.0;
    constraints.control_constraints.emplace_back("elevator", ctrl);
    EXPECT_EQ(constraints.GetNumConstraints(), 2);

    // Free control constraint should not count
    VariableConstraint free_ctrl;
    free_ctrl.type = ConstraintType::kFree;
    constraints.control_constraints.emplace_back("rudder", free_ctrl);
    EXPECT_EQ(constraints.GetNumConstraints(), 2);
}

TEST(TrimConstraintsTest, AllConstraintTypesCount) {
    TrimConstraints constraints;
    constraints.alpha.type = ConstraintType::kCL;
    constraints.beta.type = ConstraintType::kCY;
    constraints.roll_rate.type = ConstraintType::kRollMoment;
    constraints.pitch_rate.type = ConstraintType::kPitchMoment;
    constraints.yaw_rate.type = ConstraintType::kYawMoment;
    EXPECT_EQ(constraints.GetNumConstraints(), 5);
}

// ============================================================================
// TrimPresets Tests
// ============================================================================

TEST(TrimPresetsTest, LevelFlightRatesFixed) {
    TrimConstraints c = TrimPresets::LevelFlight(0.5);

    EXPECT_EQ(c.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(c.alpha.target, 0.5, kTol);
    EXPECT_EQ(c.beta.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.beta.target, 0.0, kTol);
    EXPECT_EQ(c.roll_rate.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.roll_rate.target, 0.0, kTol);
    EXPECT_EQ(c.pitch_rate.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.pitch_rate.target, 0.0, kTol);
    EXPECT_EQ(c.yaw_rate.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.yaw_rate.target, 0.0, kTol);
    EXPECT_EQ(c.GetNumConstraints(), 5);
}

TEST(TrimPresetsTest, LevelTrimmedFlight) {
    TrimConstraints c = TrimPresets::LevelTrimmedFlight(0.4, "elevator");

    // Base constraints same as LevelFlight
    EXPECT_EQ(c.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(c.alpha.target, 0.4, kTol);

    // Should have one control constraint for elevator
    ASSERT_EQ(c.control_constraints.size(), 1u);
    EXPECT_EQ(c.control_constraints[0].first, "elevator");
    EXPECT_EQ(c.control_constraints[0].second.type,
              ConstraintType::kPitchMoment);
    EXPECT_NEAR(c.control_constraints[0].second.target, 0.0, kTol);
    EXPECT_EQ(c.GetNumConstraints(), 6);
}

TEST(TrimPresetsTest, SteadyTurn) {
    Scalar bank = 30.0 * kDegreesToRadians;
    TrimConstraints c = TrimPresets::SteadyTurn(0.6, bank);

    EXPECT_EQ(c.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(c.alpha.target, 0.6, kTol);
    EXPECT_EQ(c.beta.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.beta.target, 0.0, kTol);
    EXPECT_EQ(c.roll_rate.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.roll_rate.target, 0.0, kTol);
    // Pitch and yaw rates are free for turns
    EXPECT_EQ(c.pitch_rate.type, ConstraintType::kFree);
    EXPECT_EQ(c.yaw_rate.type, ConstraintType::kFree);
    EXPECT_EQ(c.GetNumConstraints(), 3);
}

TEST(TrimPresetsTest, SteadySideslipAllFields) {
    Scalar beta_val = 5.0 * kDegreesToRadians;
    TrimConstraints c = TrimPresets::SteadySideslip(0.3, beta_val);

    EXPECT_EQ(c.alpha.type, ConstraintType::kCL);
    EXPECT_NEAR(c.alpha.target, 0.3, kTol);
    EXPECT_EQ(c.beta.type, ConstraintType::kFixed);
    EXPECT_NEAR(c.beta.target, beta_val, kTol);
    EXPECT_EQ(c.roll_rate.type, ConstraintType::kFixed);
    EXPECT_EQ(c.pitch_rate.type, ConstraintType::kFixed);
    EXPECT_EQ(c.yaw_rate.type, ConstraintType::kFixed);
    EXPECT_EQ(c.GetNumConstraints(), 5);
}

// ============================================================================
// VariableConstraint Tests
// ============================================================================

TEST(VariableConstraintTest, DefaultState) {
    VariableConstraint vc;
    EXPECT_EQ(vc.type, ConstraintType::kFree);
    EXPECT_NEAR(vc.target, 0.0, kTol);
}

}  // namespace
}  // namespace pegasus
