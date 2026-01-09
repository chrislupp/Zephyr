// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/io/mass_reader.h"

namespace vlm {
namespace {

constexpr double kTestTolerance = 1.0e-9;

// ============================================================================
// MassReader Tests
// ============================================================================

class MassReaderTest : public ::testing::Test {
 protected:
  MassReader reader_;
};

TEST_F(MassReaderTest, EmptyContent) {
  MassFileData data = reader_.ReadString("");
  EXPECT_EQ(data.components.size(), 0);
  EXPECT_NEAR(data.length_unit, 1.0, kTestTolerance);
  EXPECT_NEAR(data.mass_unit, 1.0, kTestTolerance);
  EXPECT_NEAR(data.time_unit, 1.0, kTestTolerance);
}

TEST_F(MassReaderTest, OnlyComments) {
  std::string content = R"(
# This is a comment
! Another comment
  # Indented comment
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_EQ(data.components.size(), 0);
}

TEST_F(MassReaderTest, UnitDefinitions) {
  std::string content = R"(
Lunit = 0.0254 m
Munit = 0.001 kg
Tunit = 1.0 s
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_NEAR(data.length_unit, 0.0254, kTestTolerance);
  EXPECT_NEAR(data.mass_unit, 0.001, kTestTolerance);
  EXPECT_NEAR(data.time_unit, 1.0, kTestTolerance);
}

TEST_F(MassReaderTest, ConstantDefinitions) {
  std::string content = R"(
g = 9.81
rho = 1.225
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_NEAR(data.gravity, 9.81, kTestTolerance);
  EXPECT_NEAR(data.density, 1.225, kTestTolerance);
}

TEST_F(MassReaderTest, SingleComponent) {
  std::string content = R"(
# mass   x     y     z     Ixx   Iyy   Izz
  100.0  0.5   0.0   0.0   10.0  50.0  55.0  ! fuselage
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 1);

  const MassComponent& comp = data.components[0];
  EXPECT_NEAR(comp.mass, 100.0, kTestTolerance);
  EXPECT_NEAR(comp.position.x(), 0.5, kTestTolerance);
  EXPECT_NEAR(comp.position.y(), 0.0, kTestTolerance);
  EXPECT_NEAR(comp.position.z(), 0.0, kTestTolerance);
  EXPECT_NEAR(comp.Ixx, 10.0, kTestTolerance);
  EXPECT_NEAR(comp.Iyy, 50.0, kTestTolerance);
  EXPECT_NEAR(comp.Izz, 55.0, kTestTolerance);
  EXPECT_EQ(comp.name, "fuselage");
}

TEST_F(MassReaderTest, ComponentWithProductsOfInertia) {
  std::string content = R"(
  100.0  0.5  0.0  0.0  10.0  50.0  55.0  1.0  2.0  3.0  ! with POI
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 1);

  const MassComponent& comp = data.components[0];
  EXPECT_NEAR(comp.Ixy, 1.0, kTestTolerance);
  EXPECT_NEAR(comp.Ixz, 2.0, kTestTolerance);
  EXPECT_NEAR(comp.Iyz, 3.0, kTestTolerance);
}

TEST_F(MassReaderTest, MultipleComponents) {
  std::string content = R"(
Lunit = 1.0 m
Munit = 1.0 kg

  100.0  0.5   0.0   0.0   10.0  50.0  55.0  ! fuselage
   50.0  0.3   2.0   0.1    5.0   1.0   5.0  ! right wing
   50.0  0.3  -2.0   0.1    5.0   1.0   5.0  ! left wing
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 3);
  EXPECT_EQ(data.components[0].name, "fuselage");
  EXPECT_EQ(data.components[1].name, "right wing");
  EXPECT_EQ(data.components[2].name, "left wing");
}

TEST_F(MassReaderTest, ScaleFactors) {
  std::string content = R"(
*  2.0  3.0  4.0
  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! scaled component
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 1);

  // Position should be scaled by (2, 3, 4)
  const MassComponent& comp = data.components[0];
  EXPECT_NEAR(comp.position.x(), 2.0, kTestTolerance);
  EXPECT_NEAR(comp.position.y(), 3.0, kTestTolerance);
  EXPECT_NEAR(comp.position.z(), 4.0, kTestTolerance);
}

TEST_F(MassReaderTest, OffsetFactors) {
  std::string content = R"(
+  10.0  20.0  30.0
  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! offset component
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 1);

  // Position should be offset by (10, 20, 30)
  const MassComponent& comp = data.components[0];
  EXPECT_NEAR(comp.position.x(), 11.0, kTestTolerance);
  EXPECT_NEAR(comp.position.y(), 21.0, kTestTolerance);
  EXPECT_NEAR(comp.position.z(), 31.0, kTestTolerance);
}

TEST_F(MassReaderTest, ScaleAndOffset) {
  std::string content = R"(
*  2.0  2.0  2.0
+  5.0  5.0  5.0
  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! scaled and offset
)";
  MassFileData data = reader_.ReadString(content);
  ASSERT_EQ(data.components.size(), 1);

  // Position = (1*2 + 5, 1*2 + 5, 1*2 + 5) = (7, 7, 7)
  const MassComponent& comp = data.components[0];
  EXPECT_NEAR(comp.position.x(), 7.0, kTestTolerance);
  EXPECT_NEAR(comp.position.y(), 7.0, kTestTolerance);
  EXPECT_NEAR(comp.position.z(), 7.0, kTestTolerance);
}

TEST_F(MassReaderTest, TotalMass) {
  std::string content = R"(
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp1
   50.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp2
   25.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp3
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_NEAR(data.GetTotalMass(), 175.0, kTestTolerance);
}

TEST_F(MassReaderTest, TotalMassWithUnits) {
  std::string content = R"(
Munit = 0.001 kg
  1000.0  0.0  0.0  0.0  0.0  0.0  0.0  ! 1000 grams = 1 kg
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_NEAR(data.GetTotalMass(), 1.0, kTestTolerance);
}

TEST_F(MassReaderTest, CenterOfGravity) {
  std::string content = R"(
Lunit = 1.0 m
Munit = 1.0 kg
  50.0  0.0  0.0  0.0  0.0  0.0  0.0  ! at origin
  50.0  2.0  0.0  0.0  0.0  0.0  0.0  ! at x=2
)";
  MassFileData data = reader_.ReadString(content);
  Vec3 cg = data.GetCenterOfGravity();

  // CG should be at x=1.0 (midpoint of equal masses)
  EXPECT_NEAR(cg.x(), 1.0, kTestTolerance);
  EXPECT_NEAR(cg.y(), 0.0, kTestTolerance);
  EXPECT_NEAR(cg.z(), 0.0, kTestTolerance);
}

TEST_F(MassReaderTest, CenterOfGravityWeighted) {
  std::string content = R"(
Lunit = 1.0 m
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  0.0  0.0  0.0  ! 100 kg at origin
   50.0  3.0  0.0  0.0  0.0  0.0  0.0  ! 50 kg at x=3
)";
  MassFileData data = reader_.ReadString(content);
  Vec3 cg = data.GetCenterOfGravity();

  // CG = (100*0 + 50*3) / 150 = 150/150 = 1.0
  EXPECT_NEAR(cg.x(), 1.0, kTestTolerance);
}

TEST_F(MassReaderTest, ComputeTotalProperties) {
  std::string content = R"(
Lunit = 1.0 m
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  10.0  20.0  25.0  ! single component at origin
)";
  MassFileData data = reader_.ReadString(content);
  MassProperties props = data.ComputeTotalProperties();

  EXPECT_NEAR(props.mass, 100.0, kTestTolerance);
  EXPECT_NEAR(props.cg.x(), 0.0, kTestTolerance);
  EXPECT_NEAR(props.cg.y(), 0.0, kTestTolerance);
  EXPECT_NEAR(props.cg.z(), 0.0, kTestTolerance);
  EXPECT_NEAR(props.Ixx, 10.0, kTestTolerance);
  EXPECT_NEAR(props.Iyy, 20.0, kTestTolerance);
  EXPECT_NEAR(props.Izz, 25.0, kTestTolerance);
}

TEST_F(MassReaderTest, ParallelAxisTheorem) {
  // Two point masses at y = +/- 1, each 50 kg
  // Total Ixx should include parallel axis contribution: 2 * (0 + 50 * 1^2) = 100
  std::string content = R"(
Lunit = 1.0 m
Munit = 1.0 kg
  50.0  0.0   1.0  0.0  0.0  0.0  0.0  ! at y=+1
  50.0  0.0  -1.0  0.0  0.0  0.0  0.0  ! at y=-1
)";
  MassFileData data = reader_.ReadString(content);
  MassProperties props = data.ComputeTotalProperties();

  EXPECT_NEAR(props.mass, 100.0, kTestTolerance);
  EXPECT_NEAR(props.cg.x(), 0.0, kTestTolerance);
  EXPECT_NEAR(props.cg.y(), 0.0, kTestTolerance);

  // Ixx = sum(m * (y^2 + z^2)) = 2 * 50 * (1 + 0) = 100
  EXPECT_NEAR(props.Ixx, 100.0, kTestTolerance);
  // Iyy = sum(m * (x^2 + z^2)) = 0 (both at origin in x-z)
  EXPECT_NEAR(props.Iyy, 0.0, kTestTolerance);
  // Izz = sum(m * (x^2 + y^2)) = 2 * 50 * (0 + 1) = 100
  EXPECT_NEAR(props.Izz, 100.0, kTestTolerance);
}

TEST_F(MassReaderTest, EmptyComponentsGiveZeroMass) {
  MassFileData data;
  EXPECT_NEAR(data.GetTotalMass(), 0.0, kTestTolerance);
  EXPECT_NEAR(data.GetCenterOfGravity().norm(), 0.0, kTestTolerance);
}

TEST_F(MassReaderTest, CaseInsensitiveKeywords) {
  std::string content = R"(
LUNIT = 2.0 m
munit = 3.0 kg
TuNiT = 4.0 s
G = 10.0
RHO = 1.5
)";
  MassFileData data = reader_.ReadString(content);
  EXPECT_NEAR(data.length_unit, 2.0, kTestTolerance);
  EXPECT_NEAR(data.mass_unit, 3.0, kTestTolerance);
  EXPECT_NEAR(data.time_unit, 4.0, kTestTolerance);
  EXPECT_NEAR(data.gravity, 10.0, kTestTolerance);
  EXPECT_NEAR(data.density, 1.5, kTestTolerance);
}

TEST_F(MassReaderTest, CompleteAVLStyleFile) {
  // Typical AVL mass file format
  std::string content = R"(
#-------------------------------------------------
#  Dimensional unit and calculation setup
#-------------------------------------------------
Lunit = 0.0254 m
Munit = 0.001 kg
Tunit = 1.0 s

g = 9.81
rho = 1.225

#-------------------------------------------------
#  Mass data
#-------------------------------------------------
#  mass    x       y       z     Ixx     Iyy     Izz
   500.0   10.0    0.0     0.0   100.0   200.0   250.0  ! fuselage
   100.0    5.0   20.0     1.0    10.0    50.0    55.0  ! right wing
   100.0    5.0  -20.0     1.0    10.0    50.0    55.0  ! left wing
    50.0   40.0    0.0     2.0     5.0    10.0    12.0  ! horizontal tail
    30.0   40.0    0.0     5.0     2.0     8.0     8.0  ! vertical tail
)";
  MassFileData data = reader_.ReadString(content);

  // Check units
  EXPECT_NEAR(data.length_unit, 0.0254, kTestTolerance);
  EXPECT_NEAR(data.mass_unit, 0.001, kTestTolerance);

  // Check constants
  EXPECT_NEAR(data.gravity, 9.81, kTestTolerance);
  EXPECT_NEAR(data.density, 1.225, kTestTolerance);

  // Check components
  ASSERT_EQ(data.components.size(), 5);

  // Total mass in SI: (500 + 100 + 100 + 50 + 30) * 0.001 = 0.78 kg
  EXPECT_NEAR(data.GetTotalMass(), 0.78, 1.0e-6);
}

}  // namespace
}  // namespace vlm
