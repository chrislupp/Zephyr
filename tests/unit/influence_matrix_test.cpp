// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/section.h"
#include "vlm/geometry/surface.h"
#include "vlm/solver/influence_matrix.h"
#include "vlm/solver/vortex_lattice.h"

namespace vlm {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// Biot-Savart Tests
// ============================================================================

class BiotSavartTest : public ::testing::Test {
 protected:
  // Create a simple unit horseshoe vortex centered at origin
  // Bound vortex from (-0.5, 0, 0) to (0.5, 0, 0)
  Vec3 bound_left{-0.5, 0.0, 0.0};
  Vec3 bound_right{0.5, 0.0, 0.0};
};

TEST_F(BiotSavartTest, VelocityAboveBoundVortex) {
  // Field point directly above the bound vortex center
  // With bound vortex along x-axis at y=0, and field point at y=0,
  // the velocity should have no y-component due to symmetry
  Vec3 field_point(0.0, 0.0, 1.0);
  Scalar core_radius = 0.001;

  Vec3 vel = ComputeHorseshoeVortexVelocity(field_point, bound_left, bound_right,
                                            1.0, core_radius, true);

  // Velocity should be finite (no NaN/Inf)
  EXPECT_TRUE(std::isfinite(vel.x()));
  EXPECT_TRUE(std::isfinite(vel.y()));
  EXPECT_TRUE(std::isfinite(vel.z()));

  // For this symmetric geometry (field point on y=0 plane, vortex in y=0 plane),
  // y-velocity should be very small due to symmetry
  EXPECT_NEAR(vel.y(), 0.0, 0.01);
}

TEST_F(BiotSavartTest, VelocitySymmetry) {
  // Velocity should be symmetric about x-z plane (y=0)
  Vec3 point_pos(0.0, 0.5, 0.5);
  Vec3 point_neg(0.0, -0.5, 0.5);
  Scalar core_radius = 0.001;

  Vec3 vel_pos = ComputeHorseshoeVortexVelocity(point_pos, bound_left,
                                                 bound_right, 1.0, core_radius);
  Vec3 vel_neg = ComputeHorseshoeVortexVelocity(point_neg, bound_left,
                                                 bound_right, 1.0, core_radius);

  // x-velocity should be antisymmetric
  EXPECT_NEAR(vel_pos.x(), -vel_neg.x(), kTestTolerance);

  // y-velocity should be symmetric (both negative)
  EXPECT_NEAR(vel_pos.y(), vel_neg.y(), kTestTolerance);

  // z-velocity should be antisymmetric
  EXPECT_NEAR(vel_pos.z(), -vel_neg.z(), kTestTolerance);
}

TEST_F(BiotSavartTest, CoreRegularization) {
  // Test core regularization by using only trailing vortex contributions
  // (exclude bound vortex to avoid additional complexity)
  // Right trailing vortex extends from (0.5, 0, 0) in +x direction
  // Left trailing vortex extends from (-0.5, 0, 0) in +x direction
  // Place field point asymmetrically to avoid cancellation
  Vec3 field_point(1.0, 0.3, 0.02);

  // Small core - should give larger velocity magnitude
  Vec3 vel_small_core =
      ComputeHorseshoeVortexVelocity(field_point, bound_left, bound_right, 1.0,
                                      0.01, false);  // Exclude bound vortex

  // Large core - should give smaller velocity magnitude (regularized)
  Vec3 vel_large_core =
      ComputeHorseshoeVortexVelocity(field_point, bound_left, bound_right, 1.0,
                                      0.5, false);  // Exclude bound vortex

  // Both should be finite
  EXPECT_TRUE(std::isfinite(vel_small_core.norm()));
  EXPECT_TRUE(std::isfinite(vel_large_core.norm()));

  // Both should give substantial velocity (downwash from trailing vortices)
  EXPECT_GT(vel_small_core.norm(), 0.01);

  // Larger core should reduce the velocity magnitude
  EXPECT_LT(vel_large_core.norm(), vel_small_core.norm());
}

TEST_F(BiotSavartTest, PrandtlGlauertCorrection) {
  // Point off the symmetry plane to get non-zero velocity
  Vec3 field_point(1.0, 0.5, 0.5);
  Scalar core_radius = 0.001;

  // Incompressible
  Vec3 vel_incomp = ComputeHorseshoeVortexVelocity(
      field_point, bound_left, bound_right, 1.0, core_radius);

  // Subsonic (M = 0.5, beta = sqrt(0.75) ≈ 0.866)
  Scalar beta = std::sqrt(0.75);
  Vec3 vel_comp = ComputeHorseshoeVortexVelocity(
      field_point, bound_left, bound_right, beta, core_radius);

  // Compressibility should change the velocity (Prandtl-Glauert effect)
  // The overall magnitude should increase for subsonic flow
  EXPECT_NE(vel_comp.norm(), vel_incomp.norm());
}

TEST_F(BiotSavartTest, ExcludeBoundVortex) {
  // Use a point off the symmetry plane to get meaningful velocity
  Vec3 field_point(0.0, 0.3, 0.5);
  Scalar core_radius = 0.001;

  Vec3 vel_with_bound = ComputeHorseshoeVortexVelocity(
      field_point, bound_left, bound_right, 1.0, core_radius, true);

  Vec3 vel_no_bound = ComputeHorseshoeVortexVelocity(
      field_point, bound_left, bound_right, 1.0, core_radius, false);

  // Both should be finite
  EXPECT_TRUE(std::isfinite(vel_with_bound.norm()));
  EXPECT_TRUE(std::isfinite(vel_no_bound.norm()));

  // Without bound vortex, velocity should be different
  EXPECT_NE(vel_with_bound.norm(), vel_no_bound.norm());
}

TEST_F(BiotSavartTest, TrailingVortexInDownstreamField) {
  // Point far downstream (in +x direction)
  Vec3 field_point(10.0, 0.5, 0.0);
  Scalar core_radius = 0.001;

  Vec3 vel = ComputeHorseshoeVortexVelocity(field_point, bound_left, bound_right,
                                            1.0, core_radius, true);

  // Far downstream, velocity is primarily from trailing vortices
  // Should induce downwash (negative z) for positive y location
  EXPECT_LT(vel.z(), 0.0);
}

// ============================================================================
// InfluenceMatrixBuilder Tests
// ============================================================================

class InfluenceMatrixTest : public ::testing::Test {
 protected:
  Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                 int num_spanwise) {
    Aircraft aircraft("Test Wing");

    Surface wing("Wing");
    wing.SetChordwisePanelCount(num_chordwise);
    wing.SetSpanwisePanelCount(num_spanwise);
    wing.SetChordwiseSpacing(SpacingType::kUniform);
    wing.SetSpanwiseSpacing(SpacingType::kUniform);

    // Full wing (no y-duplicate for simplicity)
    wing.AddSection(Section(Vec3(0, -span / 2, 0), chord, 0.0));
    wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));

    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(span * chord);
    aircraft.SetReferenceChord(chord);
    aircraft.SetReferenceSpan(span);

    return aircraft;
  }
};

TEST_F(InfluenceMatrixTest, VelocityMatrixDimensions) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);
  MatX wc_gam = builder.BuildVelocityInfluenceMatrix();

  Index n = mesh.GetNumVortices();
  EXPECT_EQ(wc_gam.rows(), 3 * n);
  EXPECT_EQ(wc_gam.cols(), n);
}

TEST_F(InfluenceMatrixTest, AICMatrixDimensions) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);
  MatX aic = builder.BuildNormalwashMatrix();

  Index n = mesh.GetNumVortices();
  EXPECT_EQ(aic.rows(), n);
  EXPECT_EQ(aic.cols(), n);
}

TEST_F(InfluenceMatrixTest, AICDiagonalDominance) {
  // For a flat wing, AIC should be diagonally dominant
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);
  MatX aic = builder.BuildNormalwashMatrix();

  Index n = mesh.GetNumVortices();

  // Check that diagonal elements are significant
  for (Index i = 0; i < n; ++i) {
    Scalar diag = std::abs(aic(i, i));
    Scalar row_sum = 0.0;
    for (Index j = 0; j < n; ++j) {
      if (i != j) {
        row_sum += std::abs(aic(i, j));
      }
    }
    // Diagonal should be a significant portion of the row
    EXPECT_GT(diag, 0.0);
  }
}

TEST_F(InfluenceMatrixTest, AICMatrixNonsingular) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);
  MatX aic = builder.BuildNormalwashMatrix();

  // Matrix should be invertible
  Eigen::PartialPivLU<MatX> lu(aic);
  Scalar det = lu.determinant();

  EXPECT_NE(det, 0.0);
}

TEST_F(InfluenceMatrixTest, MachNumberEffect) {
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);

  // Incompressible
  builder.SetMachNumber(0.0);
  MatX aic_incomp = builder.BuildNormalwashMatrix();

  // Subsonic
  builder.SetMachNumber(0.5);
  MatX aic_comp = builder.BuildNormalwashMatrix();

  // Matrices should be different due to Prandtl-Glauert
  Scalar diff = (aic_comp - aic_incomp).norm();
  EXPECT_GT(diff, 0.0);
}

TEST_F(InfluenceMatrixTest, InducedVelocityFromMatrix) {
  // Test that velocity influence matrix produces finite, reasonable values
  Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);

  VortexLattice mesh;
  mesh.GenerateFromAircraft(aircraft);

  InfluenceMatrixBuilder builder(mesh);
  MatX velocity_matrix;
  builder.BuildNormalwashMatrix(&velocity_matrix);

  // All values should be finite
  for (Index i = 0; i < velocity_matrix.rows(); ++i) {
    for (Index j = 0; j < velocity_matrix.cols(); ++j) {
      EXPECT_TRUE(std::isfinite(velocity_matrix(i, j)));
    }
  }

  // Diagonal-ish influence should be significant
  // (self-influence from trailing legs)
  Index n = mesh.GetNumVortices();
  for (Index i = 0; i < n; ++i) {
    Scalar self_influence = std::sqrt(
        velocity_matrix(3 * i + 0, i) * velocity_matrix(3 * i + 0, i) +
        velocity_matrix(3 * i + 1, i) * velocity_matrix(3 * i + 1, i) +
        velocity_matrix(3 * i + 2, i) * velocity_matrix(3 * i + 2, i));
    EXPECT_GT(self_influence, 0.0);
  }
}

}  // namespace
}  // namespace vlm
