// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <gtest/gtest.h>

#include <cmath>

#include "vlm/ad/enzyme_geometry.h"
#include "vlm/core/types.h"

namespace vlm {
namespace ad {
namespace {

// ============================================================================
// Test Fixture
// ============================================================================

class EnzymeGeometryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Simple rectangular wing section data
    // Left section at y=0, right section at y=2
    le_left_[0] = 0.0;  // x
    le_left_[1] = 0.0;  // y
    le_left_[2] = 0.0;  // z

    le_right_[0] = 0.0;  // x
    le_right_[1] = 2.0;  // y
    le_right_[2] = 0.0;  // z

    chord_left_ = 1.0;
    chord_right_ = 1.0;
    twist_left_ = 0.0;
    twist_right_ = 0.0;

    // Panel at midspan
    frac_left_ = 0.4;
    frac_center_ = 0.5;
    frac_right_ = 0.6;

    x_vortex_ = 0.25;
    x_control_ = 0.75;
  }

  double le_left_[3];
  double le_right_[3];
  double chord_left_, chord_right_;
  double twist_left_, twist_right_;
  double frac_left_, frac_center_, frac_right_;
  double x_vortex_, x_control_;
};

// ============================================================================
// Basic Kernel Tests
// ============================================================================

TEST_F(EnzymeGeometryTest, InterpolatedLEIsCorrect) {
  double le_out[3];
  compute_interpolated_le(le_left_, le_right_, 0.5, le_out);

  EXPECT_DOUBLE_EQ(le_out[0], 0.0);
  EXPECT_DOUBLE_EQ(le_out[1], 1.0);  // Midpoint
  EXPECT_DOUBLE_EQ(le_out[2], 0.0);
}

TEST_F(EnzymeGeometryTest, InterpolatedChordIsCorrect) {
  double chord = compute_interpolated_chord(1.0, 0.5, 0.5);
  EXPECT_DOUBLE_EQ(chord, 0.75);
}

TEST_F(EnzymeGeometryTest, InterpolatedTwistZeroForZeroTwist) {
  double twist = compute_interpolated_twist(0.0, 0.0, 1.0, 1.0, 0.5);
  EXPECT_NEAR(twist, 0.0, 1e-10);
}

TEST_F(EnzymeGeometryTest, InterpolatedTwistIsCorrect) {
  double twist = compute_interpolated_twist(
      10.0 * kDegreesToRadians, 10.0 * kDegreesToRadians,
      1.0, 1.0, 0.5);
  EXPECT_NEAR(twist, 10.0 * kDegreesToRadians, 1e-10);
}

TEST_F(EnzymeGeometryTest, TrefftzNormalForXYPlaneWing) {
  // Wing in XY plane: Trefftz normal should be (0, 0, 1)
  double trefftz_y, trefftz_z;
  compute_trefftz_normal(le_left_, le_right_, &trefftz_y, &trefftz_z);

  EXPECT_NEAR(trefftz_y, 0.0, 1e-10);
  EXPECT_NEAR(trefftz_z, 1.0, 1e-10);
}

TEST_F(EnzymeGeometryTest, TrefftzNormalForDihedralWing) {
  // Wing with 45 deg dihedral: LE goes from (0,0,0) to (0,1,1)
  double le_right_dihedral[3] = {0.0, 1.0, 1.0};
  double trefftz_y, trefftz_z;
  compute_trefftz_normal(le_left_, le_right_dihedral, &trefftz_y, &trefftz_z);

  // Expected: perpendicular to (0,1,1) in YZ plane
  // That's proportional to (-1, 1) normalized = (-1/sqrt(2), 1/sqrt(2))
  double expected_y = -1.0 / std::sqrt(2.0);
  double expected_z = 1.0 / std::sqrt(2.0);

  EXPECT_NEAR(trefftz_y, expected_y, 1e-10);
  EXPECT_NEAR(trefftz_z, expected_z, 1e-10);
}

TEST_F(EnzymeGeometryTest, ControlPointPosition) {
  double le[3] = {1.0, 2.0, 3.0};
  double cp[3];
  compute_control_point(le, 1.0, 0.75, cp);

  EXPECT_DOUBLE_EQ(cp[0], 1.75);  // x + 0.75*chord
  EXPECT_DOUBLE_EQ(cp[1], 2.0);   // y unchanged
  EXPECT_DOUBLE_EQ(cp[2], 3.0);   // z unchanged
}

TEST_F(EnzymeGeometryTest, PanelNormalForZeroTwist) {
  // For zero twist, chord direction is (1, 0, 0)
  // Bound direction along y is (0, 1, 0)
  // Normal = (1,0,0) x (0,1,0) = (0, 0, 1)
  double bound_dir[3] = {0.0, 1.0, 0.0};
  double normal[3];
  compute_panel_normal(0.0, 0.0, 1.0, bound_dir, normal);

  EXPECT_NEAR(normal[0], 0.0, 1e-10);
  EXPECT_NEAR(normal[1], 0.0, 1e-10);
  EXPECT_NEAR(normal[2], 1.0, 1e-10);
}

TEST_F(EnzymeGeometryTest, PanelNormalWithTwist) {
  // 5 deg twist about (0, 0, 1) Trefftz normal
  double twist = 5.0 * kDegreesToRadians;
  double bound_dir[3] = {0.0, 1.0, 0.0};
  double normal[3];
  compute_panel_normal(twist, 0.0, 1.0, bound_dir, normal);

  // Normal should tilt aft (positive x component) with nose-up twist
  EXPECT_GT(normal[0], 0.0);  // Tilted aft
  EXPECT_NEAR(normal[1], 0.0, 1e-10);  // No sideways tilt
  EXPECT_GT(normal[2], 0.0);  // Still mostly up
}

// ============================================================================
// Full Panel Geometry Computation
// ============================================================================

TEST_F(EnzymeGeometryTest, FullPanelGeometry) {
  double bv_left[3], bv_right[3], cp[3], normal[3];

  compute_panel_geometry(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      bv_left, bv_right, cp, normal);

  // Check bound vortex left position
  // LE at frac=0.4 is (0, 0.8, 0), chord=1.0, vortex at x=0.25
  EXPECT_NEAR(bv_left[0], 0.25, 1e-10);
  EXPECT_NEAR(bv_left[1], 0.8, 1e-10);
  EXPECT_NEAR(bv_left[2], 0.0, 1e-10);

  // Check control point
  // LE at frac=0.5 is (0, 1.0, 0), chord=1.0, control at x=0.75
  EXPECT_NEAR(cp[0], 0.75, 1e-10);
  EXPECT_NEAR(cp[1], 1.0, 1e-10);
  EXPECT_NEAR(cp[2], 0.0, 1e-10);

  // Normal should be (0, 0, 1) for zero twist XY plane wing
  EXPECT_NEAR(normal[0], 0.0, 1e-10);
  EXPECT_NEAR(normal[1], 0.0, 1e-10);
  EXPECT_NEAR(normal[2], 1.0, 1e-10);
}

// ============================================================================
// Finite Difference Gradient Tests
// ============================================================================

TEST_F(EnzymeGeometryTest, ChordGradientFD) {
  double d_bv_left[3], d_bv_right[3], d_cp[3], d_normal[3];

  fd_panel_geometry_chord_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,  // left section
      d_bv_left, d_bv_right, d_cp, d_normal);

  // d(bv_left.x)/d(chord_left) should be positive (vortex moves aft with chord)
  // For left section at frac=0.4, contribution is (1-0.4)*x_vortex = 0.6*0.25 = 0.15
  EXPECT_NEAR(d_bv_left[0], 0.6 * x_vortex_, 1e-5);

  // d(cp.x)/d(chord_left) at center (frac=0.5) is 0.5*x_control = 0.375
  EXPECT_NEAR(d_cp[0], 0.5 * x_control_, 1e-5);

  // y and z positions don't depend on chord
  EXPECT_NEAR(d_bv_left[1], 0.0, 1e-8);
  EXPECT_NEAR(d_bv_left[2], 0.0, 1e-8);
}

TEST_F(EnzymeGeometryTest, TwistGradientFD) {
  double d_bv_left[3], d_bv_right[3], d_cp[3], d_normal[3];

  fd_panel_geometry_twist_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,  // left section
      d_bv_left, d_bv_right, d_cp, d_normal);

  // Positions don't depend on twist
  EXPECT_NEAR(d_bv_left[0], 0.0, 1e-8);
  EXPECT_NEAR(d_bv_left[1], 0.0, 1e-8);
  EXPECT_NEAR(d_bv_left[2], 0.0, 1e-8);

  // Normal DOES depend on twist
  // For small twist, d(normal.x)/d(twist) should be positive (tilts aft)
  EXPECT_GT(d_normal[0], 0.0);
}

TEST_F(EnzymeGeometryTest, LEPositionGradientFD) {
  double d_bv_left[3], d_bv_right[3], d_cp[3], d_normal[3];

  // Gradient w.r.t. x-position of left section
  fd_panel_geometry_le_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,  // left section
      0,  // x component
      d_bv_left, d_bv_right, d_cp, d_normal);

  // d(bv_left.x)/d(le_left.x) at frac=0.4 is 1-0.4 = 0.6
  EXPECT_NEAR(d_bv_left[0], 0.6, 1e-5);

  // d(cp.x)/d(le_left.x) at frac=0.5 is 0.5
  EXPECT_NEAR(d_cp[0], 0.5, 1e-5);
}

// ============================================================================
// Enzyme Gradient Tests (conditional on VLM_USE_ENZYME)
// ============================================================================

#ifdef VLM_USE_ENZYME

TEST_F(EnzymeGeometryTest, EnzymeChordGradientMatchesFD) {
  double d_bv_left_enzyme[3], d_bv_right_enzyme[3], d_cp_enzyme[3], d_normal_enzyme[3];
  double d_bv_left_fd[3], d_bv_right_fd[3], d_cp_fd[3], d_normal_fd[3];

  enzyme_panel_geometry_chord_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,  // left section
      d_bv_left_enzyme, d_bv_right_enzyme, d_cp_enzyme, d_normal_enzyme);

  fd_panel_geometry_chord_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,  // left section
      d_bv_left_fd, d_bv_right_fd, d_cp_fd, d_normal_fd);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(d_bv_left_enzyme[i], d_bv_left_fd[i], 1e-6)
        << "Mismatch in d_bv_left[" << i << "]";
    EXPECT_NEAR(d_bv_right_enzyme[i], d_bv_right_fd[i], 1e-6)
        << "Mismatch in d_bv_right[" << i << "]";
    EXPECT_NEAR(d_cp_enzyme[i], d_cp_fd[i], 1e-6)
        << "Mismatch in d_cp[" << i << "]";
    EXPECT_NEAR(d_normal_enzyme[i], d_normal_fd[i], 1e-6)
        << "Mismatch in d_normal[" << i << "]";
  }
}

TEST_F(EnzymeGeometryTest, EnzymeTwistGradientMatchesFD) {
  double d_bv_left_enzyme[3], d_bv_right_enzyme[3], d_cp_enzyme[3], d_normal_enzyme[3];
  double d_bv_left_fd[3], d_bv_right_fd[3], d_cp_fd[3], d_normal_fd[3];

  enzyme_panel_geometry_twist_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,
      d_bv_left_enzyme, d_bv_right_enzyme, d_cp_enzyme, d_normal_enzyme);

  fd_panel_geometry_twist_gradient(
      le_left_, le_right_,
      chord_left_, chord_right_,
      twist_left_, twist_right_,
      frac_left_, frac_center_, frac_right_,
      x_vortex_, x_control_,
      0,
      d_bv_left_fd, d_bv_right_fd, d_cp_fd, d_normal_fd);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(d_bv_left_enzyme[i], d_bv_left_fd[i], 1e-6);
    EXPECT_NEAR(d_bv_right_enzyme[i], d_bv_right_fd[i], 1e-6);
    EXPECT_NEAR(d_cp_enzyme[i], d_cp_fd[i], 1e-6);
    EXPECT_NEAR(d_normal_enzyme[i], d_normal_fd[i], 1e-6);
  }
}

TEST_F(EnzymeGeometryTest, EnzymeLEGradientMatchesFD) {
  double d_bv_left_enzyme[3], d_bv_right_enzyme[3], d_cp_enzyme[3], d_normal_enzyme[3];
  double d_bv_left_fd[3], d_bv_right_fd[3], d_cp_fd[3], d_normal_fd[3];

  for (int component = 0; component < 3; ++component) {
    enzyme_panel_geometry_le_gradient(
        le_left_, le_right_,
        chord_left_, chord_right_,
        twist_left_, twist_right_,
        frac_left_, frac_center_, frac_right_,
        x_vortex_, x_control_,
        0, component,
        d_bv_left_enzyme, d_bv_right_enzyme, d_cp_enzyme, d_normal_enzyme);

    fd_panel_geometry_le_gradient(
        le_left_, le_right_,
        chord_left_, chord_right_,
        twist_left_, twist_right_,
        frac_left_, frac_center_, frac_right_,
        x_vortex_, x_control_,
        0, component,
        d_bv_left_fd, d_bv_right_fd, d_cp_fd, d_normal_fd);

    for (int i = 0; i < 3; ++i) {
      EXPECT_NEAR(d_bv_left_enzyme[i], d_bv_left_fd[i], 1e-6)
          << "component=" << component << " i=" << i;
      EXPECT_NEAR(d_cp_enzyme[i], d_cp_fd[i], 1e-6)
          << "component=" << component << " i=" << i;
    }
  }
}

#endif  // VLM_USE_ENZYME

}  // namespace
}  // namespace ad
}  // namespace vlm
