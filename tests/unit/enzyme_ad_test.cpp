// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/ad/enzyme_interface.h"
#include "vlm/ad/enzyme_forces.h"
#include "vlm/ad/enzyme_biot_savart.h"
#include "vlm/ad/enzyme_mesh.h"
#include "vlm/ad/enzyme_trefftz.h"

namespace vlm {
namespace {

#ifdef VLM_USE_ENZYME
using namespace ad;  // Import Enzyme symbols (__enzyme_autodiff, enzyme_dup, etc.)
#endif

// ============================================================================
// Test Functions for Differentiation
// ============================================================================

// Simple scalar function: f(x) = x^2
// df/dx = 2x
double square(double x) { return x * x; }

// Two-argument function: f(x,y) = x^2 + y^2
// df/dx = 2x, df/dy = 2y
double sum_of_squares(double x, double y) { return x * x + y * y; }

// Product function: f(x,y) = x * y
// df/dx = y, df/dy = x
double product(double x, double y) { return x * y; }

// Exponential: f(x) = exp(x)
// df/dx = exp(x)
double exponential(double x) { return std::exp(x); }

// Trigonometric: f(x) = sin(x)
// df/dx = cos(x)
double sine(double x) { return std::sin(x); }

// Composite: f(x) = sin(x^2)
// df/dx = 2x * cos(x^2)
double sin_of_square(double x) { return std::sin(x * x); }

// Rosenbrock function (2D): f(x,y) = (1-x)^2 + 100*(y-x^2)^2
// df/dx = -2*(1-x) - 400*x*(y-x^2)
// df/dy = 200*(y-x^2)
double rosenbrock(double x, double y) {
  double a = 1.0 - x;
  double b = y - x * x;
  return a * a + 100.0 * b * b;
}

// Array function: computes sum of squared elements
// d/dx_i = 2*x_i
double sum_squared_array(const double* x, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += x[i] * x[i];
  }
  return sum;
}

// ============================================================================
// Basic Enzyme Tests
// ============================================================================

class EnzymeBasicTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-10;
};

TEST_F(EnzymeBasicTest, IsAvailableReturnsCorrectValue) {
#ifdef VLM_USE_ENZYME
  EXPECT_TRUE(ad::IsEnzymeAvailable());
#else
  EXPECT_FALSE(ad::IsEnzymeAvailable());
#endif
}

#ifdef VLM_USE_ENZYME

TEST_F(EnzymeBasicTest, SquareGradient) {
  // f(x) = x^2, df/dx = 2x
  // At x=3: df/dx = 6
  double x = 3.0;
  double grad = __enzyme_autodiff<double>((void*)square, x);

  EXPECT_NEAR(grad, 6.0, kTolerance);
}

TEST_F(EnzymeBasicTest, SquareGradientAtZero) {
  // At x=0: df/dx = 0
  double x = 0.0;
  double grad = __enzyme_autodiff<double>((void*)square, x);

  EXPECT_NEAR(grad, 0.0, kTolerance);
}

TEST_F(EnzymeBasicTest, SquareGradientNegative) {
  // At x=-2: df/dx = -4
  double x = -2.0;
  double grad = __enzyme_autodiff<double>((void*)square, x);

  EXPECT_NEAR(grad, -4.0, kTolerance);
}

TEST_F(EnzymeBasicTest, ExponentialGradient) {
  // f(x) = exp(x), df/dx = exp(x)
  // At x=1: df/dx = e
  double x = 1.0;
  double grad = __enzyme_autodiff<double>((void*)exponential, x);

  EXPECT_NEAR(grad, std::exp(1.0), kTolerance);
}

TEST_F(EnzymeBasicTest, SineGradient) {
  // f(x) = sin(x), df/dx = cos(x)
  // At x=pi/4: df/dx = cos(pi/4) = sqrt(2)/2
  double x = M_PI / 4.0;
  double grad = __enzyme_autodiff<double>((void*)sine, x);

  EXPECT_NEAR(grad, std::cos(M_PI / 4.0), kTolerance);
}

TEST_F(EnzymeBasicTest, CompositeGradient) {
  // f(x) = sin(x^2), df/dx = 2x*cos(x^2)
  // At x=1: df/dx = 2*cos(1)
  double x = 1.0;
  double grad = __enzyme_autodiff<double>((void*)sin_of_square, x);

  EXPECT_NEAR(grad, 2.0 * std::cos(1.0), kTolerance);
}

// ============================================================================
// Multi-Variable Tests
// ============================================================================

class EnzymeMultiVariableTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-10;
};

TEST_F(EnzymeMultiVariableTest, SumOfSquaresGradient) {
  // f(x,y) = x^2 + y^2
  // At (3,4): df/dx = 6, df/dy = 8
  double x = 3.0, y = 4.0;

  // Use forward mode to compute partial derivatives
  double dx = __enzyme_fwddiff<double>((void*)sum_of_squares,
                                        enzyme_dup, x, 1.0,
                                        enzyme_dup, y, 0.0);
  double dy = __enzyme_fwddiff<double>((void*)sum_of_squares,
                                        enzyme_dup, x, 0.0,
                                        enzyme_dup, y, 1.0);

  EXPECT_NEAR(dx, 6.0, kTolerance);
  EXPECT_NEAR(dy, 8.0, kTolerance);
}

TEST_F(EnzymeMultiVariableTest, ProductGradient) {
  // f(x,y) = x*y
  // df/dx = y, df/dy = x
  // At (3,5): df/dx = 5, df/dy = 3
  double x = 3.0, y = 5.0;

  double dx = __enzyme_fwddiff<double>((void*)product,
                                        enzyme_dup, x, 1.0,
                                        enzyme_dup, y, 0.0);
  double dy = __enzyme_fwddiff<double>((void*)product,
                                        enzyme_dup, x, 0.0,
                                        enzyme_dup, y, 1.0);

  EXPECT_NEAR(dx, 5.0, kTolerance);
  EXPECT_NEAR(dy, 3.0, kTolerance);
}

TEST_F(EnzymeMultiVariableTest, RosenbrockGradient) {
  // Rosenbrock: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
  // At (1,1) which is the minimum: df/dx = 0, df/dy = 0
  // Use forward mode to compute partial derivatives
  double x = 1.0, y = 1.0;

  // Compute df/dx using forward mode (seed x with 1, y with 0)
  double dx = __enzyme_fwddiff<double>((void*)rosenbrock,
                                        enzyme_dup, x, 1.0,
                                        enzyme_dup, y, 0.0);

  double dy = __enzyme_fwddiff<double>((void*)rosenbrock,
                                        enzyme_dup, x, 0.0,
                                        enzyme_dup, y, 1.0);

  EXPECT_NEAR(dx, 0.0, kTolerance);
  EXPECT_NEAR(dy, 0.0, kTolerance);
}

TEST_F(EnzymeMultiVariableTest, RosenbrockGradientAwayFromMinimum) {
  // At (0,0): df/dx = -2, df/dy = 0
  double x = 0.0, y = 0.0;

  double dx = __enzyme_fwddiff<double>((void*)rosenbrock,
                                        enzyme_dup, x, 1.0,
                                        enzyme_dup, y, 0.0);

  double dy = __enzyme_fwddiff<double>((void*)rosenbrock,
                                        enzyme_dup, x, 0.0,
                                        enzyme_dup, y, 1.0);

  EXPECT_NEAR(dx, -2.0, kTolerance);
  EXPECT_NEAR(dy, 0.0, kTolerance);
}

TEST_F(EnzymeMultiVariableTest, ConstantArgumentNotDifferentiated) {
  // Only differentiate with respect to x, treat y as constant
  double x = 3.0, y = 5.0;

  // f(x,y) = x*y, df/dx = y = 5
  double dx = __enzyme_fwddiff<double>((void*)product,
                                        enzyme_dup, x, 1.0,
                                        enzyme_const, y);

  EXPECT_NEAR(dx, 5.0, kTolerance);
}

// ============================================================================
// Forward Mode Tests
// ============================================================================

class EnzymeForwardModeTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-10;
};

TEST_F(EnzymeForwardModeTest, SquareForwardDerivative) {
  // Forward mode: f(x) = x^2, df/dx = 2x
  // At x=3 with dx=1: returns df = 6
  double x = 3.0;
  double dx = 1.0;  // Seed tangent

  double df = __enzyme_fwddiff<double>((void*)square,
                                        enzyme_dup, x, dx);

  EXPECT_NEAR(df, 6.0, kTolerance);
}

TEST_F(EnzymeForwardModeTest, ForwardModeScaling) {
  // Forward mode with scaled tangent
  // At x=3 with dx=2: returns df = 2*6 = 12
  double x = 3.0;
  double dx = 2.0;

  double df = __enzyme_fwddiff<double>((void*)square,
                                        enzyme_dup, x, dx);

  EXPECT_NEAR(df, 12.0, kTolerance);
}

// ============================================================================
// Finite Difference Verification
// ============================================================================

class EnzymeFiniteDiffTest : public ::testing::Test {
 protected:
  static constexpr double kFDStep = 1e-7;
  static constexpr double kFDTolerance = 1e-5;

  // Compute finite difference gradient
  static double FiniteDiff(double (*f)(double), double x) {
    return (f(x + kFDStep) - f(x - kFDStep)) / (2.0 * kFDStep);
  }
};

TEST_F(EnzymeFiniteDiffTest, SquareMatchesFiniteDiff) {
  double x = 2.5;
  double enzyme_grad = __enzyme_autodiff<double>((void*)square, x);
  double fd_grad = FiniteDiff(square, x);

  EXPECT_NEAR(enzyme_grad, fd_grad, kFDTolerance);
}

TEST_F(EnzymeFiniteDiffTest, ExponentialMatchesFiniteDiff) {
  double x = 1.5;
  double enzyme_grad = __enzyme_autodiff<double>((void*)exponential, x);
  double fd_grad = FiniteDiff(exponential, x);

  EXPECT_NEAR(enzyme_grad, fd_grad, kFDTolerance);
}

TEST_F(EnzymeFiniteDiffTest, SineMatchesFiniteDiff) {
  double x = 0.7;
  double enzyme_grad = __enzyme_autodiff<double>((void*)sine, x);
  double fd_grad = FiniteDiff(sine, x);

  EXPECT_NEAR(enzyme_grad, fd_grad, kFDTolerance);
}

TEST_F(EnzymeFiniteDiffTest, CompositeMatchesFiniteDiff) {
  double x = 0.8;
  double enzyme_grad = __enzyme_autodiff<double>((void*)sin_of_square, x);
  double fd_grad = FiniteDiff(sin_of_square, x);

  EXPECT_NEAR(enzyme_grad, fd_grad, kFDTolerance);
}

// ============================================================================
// Force Kernel Tests
// ============================================================================

class EnzymeForceKernelTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-6;

  void SetUp() override {
    // Simple 2-vortex setup for testing
    num_vortices_ = 2;
    reference_area_ = 1.0;

    // Circulation values
    circulation_[0] = 0.5;
    circulation_[1] = 0.3;

    // Bound vortex endpoints (two panels side by side)
    // Panel 0: left at (-0.5, 0, 0), right at (0, 0, 0)
    bound_left_[0] = -0.5; bound_left_[1] = 0.0; bound_left_[2] = 0.0;
    bound_right_[0] = 0.0; bound_right_[1] = 0.0; bound_right_[2] = 0.0;
    // Panel 1: left at (0, 0, 0), right at (0.5, 0, 0)
    bound_left_[3] = 0.0; bound_left_[4] = 0.0; bound_left_[5] = 0.0;
    bound_right_[3] = 0.5; bound_right_[4] = 0.0; bound_right_[5] = 0.0;

    // Control points at panel centers
    control_points_[0] = -0.25; control_points_[1] = 0.0; control_points_[2] = 0.0;
    control_points_[3] = 0.25; control_points_[4] = 0.0; control_points_[5] = 0.0;

    // Moment arms (relative to origin)
    moment_arm_[0] = -0.25; moment_arm_[1] = 0.0; moment_arm_[2] = 0.0;
    moment_arm_[3] = 0.25; moment_arm_[4] = 0.0; moment_arm_[5] = 0.0;

    // Small induced velocities
    for (int i = 0; i < 6; ++i) {
      induced_velocity_[i] = 0.01 * (i + 1);
    }
  }

  int num_vortices_;
  double reference_area_;
  double circulation_[2];
  double bound_left_[6];
  double bound_right_[6];
  double control_points_[6];
  double moment_arm_[6];
  double induced_velocity_[6];
};

TEST_F(EnzymeForceKernelTest, FreestreamDirectionAtZeroAlphaBeta) {
  double v_inf[3];
  ad::compute_freestream_direction(0.0, 0.0, v_inf);

  // At alpha=0, beta=0: freestream is (1, 0, 0)
  EXPECT_NEAR(v_inf[0], 1.0, kTolerance);
  EXPECT_NEAR(v_inf[1], 0.0, kTolerance);
  EXPECT_NEAR(v_inf[2], 0.0, kTolerance);
}

TEST_F(EnzymeForceKernelTest, FreestreamDirectionAtAlpha) {
  double alpha = 0.1;  // ~5.7 degrees
  double v_inf[3];
  ad::compute_freestream_direction(alpha, 0.0, v_inf);

  // At alpha, beta=0: v_inf = (cos(alpha), 0, sin(alpha))
  EXPECT_NEAR(v_inf[0], std::cos(alpha), kTolerance);
  EXPECT_NEAR(v_inf[1], 0.0, kTolerance);
  EXPECT_NEAR(v_inf[2], std::sin(alpha), kTolerance);
}

TEST_F(EnzymeForceKernelTest, FreestreamAlphaGradientMatchesFiniteDiff) {
  double alpha = 0.1;
  double beta = 0.05;
  double h = 1e-7;

  // Analytical gradient
  double dv_dalpha[3];
  ad::compute_freestream_alpha_gradient(alpha, beta, dv_dalpha);

  // Finite difference
  double v_plus[3], v_minus[3];
  ad::compute_freestream_direction(alpha + h, beta, v_plus);
  ad::compute_freestream_direction(alpha - h, beta, v_minus);

  for (int i = 0; i < 3; ++i) {
    double fd = (v_plus[i] - v_minus[i]) / (2 * h);
    EXPECT_NEAR(dv_dalpha[i], fd, kTolerance)
        << "Mismatch at component " << i;
  }
}

TEST_F(EnzymeForceKernelTest, FreestreamBetaGradientMatchesFiniteDiff) {
  double alpha = 0.1;
  double beta = 0.05;
  double h = 1e-7;

  // Analytical gradient
  double dv_dbeta[3];
  ad::compute_freestream_beta_gradient(alpha, beta, dv_dbeta);

  // Finite difference
  double v_plus[3], v_minus[3];
  ad::compute_freestream_direction(alpha, beta + h, v_plus);
  ad::compute_freestream_direction(alpha, beta - h, v_minus);

  for (int i = 0; i < 3; ++i) {
    double fd = (v_plus[i] - v_minus[i]) / (2 * h);
    EXPECT_NEAR(dv_dbeta[i], fd, kTolerance)
        << "Mismatch at component " << i;
  }
}

TEST_F(EnzymeForceKernelTest, FreestreamAlphaGradientAtZero) {
  // At alpha=0, beta=0: ∂V/∂α = (-sin(0)*cos(0), 0, cos(0)*cos(0)) = (0, 0, 1)
  double dv_dalpha[3];
  ad::compute_freestream_alpha_gradient(0.0, 0.0, dv_dalpha);

  EXPECT_NEAR(dv_dalpha[0], 0.0, kTolerance);
  EXPECT_NEAR(dv_dalpha[1], 0.0, kTolerance);
  EXPECT_NEAR(dv_dalpha[2], 1.0, kTolerance);
}

TEST_F(EnzymeForceKernelTest, FreestreamBetaGradientAtZero) {
  // At alpha=0, beta=0: ∂V/∂β = (-cos(0)*sin(0), cos(0), -sin(0)*sin(0)) = (0, 1, 0)
  double dv_dbeta[3];
  ad::compute_freestream_beta_gradient(0.0, 0.0, dv_dbeta);

  EXPECT_NEAR(dv_dbeta[0], 0.0, kTolerance);
  EXPECT_NEAR(dv_dbeta[1], 1.0, kTolerance);
  EXPECT_NEAR(dv_dbeta[2], 0.0, kTolerance);
}

TEST_F(EnzymeForceKernelTest, BodyToStabilityAtZeroAlpha) {
  double cf_x = 0.1, cf_z = 1.0;
  double CL, CD;
  ad::body_to_stability(cf_x, cf_z, 0.0, &CL, &CD);

  // At alpha=0: CL = cf_z, CD = cf_x
  EXPECT_NEAR(CL, cf_z, kTolerance);
  EXPECT_NEAR(CD, cf_x, kTolerance);
}

TEST_F(EnzymeForceKernelTest, ForceKernelProducesNonzeroForces) {
  double alpha = 0.05;
  double beta = 0.0;
  double result[6];

  ad::compute_force_coefficients_kernel(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, result);

  // Should produce non-zero lift
  EXPECT_NE(result[0], 0.0);  // CL
}

TEST_F(EnzymeForceKernelTest, AlphaGradientMatchesFiniteDiff) {
  double alpha = 0.05;
  double beta = 0.0;

  // Enzyme gradient
  double d_result_enzyme[6];
  ad::enzyme_force_alpha_gradient(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, d_result_enzyme);

  // Finite difference gradient
  double d_result_fd[6];
  ad::fd_force_alpha_gradient(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, d_result_fd);

  // Compare
  for (int i = 0; i < 6; ++i) {
    EXPECT_NEAR(d_result_enzyme[i], d_result_fd[i], kTolerance)
        << "Mismatch at output index " << i;
  }
}

TEST_F(EnzymeForceKernelTest, BetaGradientMatchesFiniteDiff) {
  double alpha = 0.05;
  double beta = 0.02;

  // Enzyme gradient
  double d_result_enzyme[6];
  ad::enzyme_force_beta_gradient(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, d_result_enzyme);

  // Finite difference gradient
  double d_result_fd[6];
  ad::fd_force_beta_gradient(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, d_result_fd);

  // Compare
  for (int i = 0; i < 6; ++i) {
    EXPECT_NEAR(d_result_enzyme[i], d_result_fd[i], kTolerance)
        << "Mismatch at output index " << i;
  }
}

TEST_F(EnzymeForceKernelTest, CirculationGradientMatchesFiniteDiff) {
  double alpha = 0.05;
  double beta = 0.0;
  double h = 1e-7;

  // Test gradient of CL (output_index=0) w.r.t. circulation
  double d_circulation_enzyme[2];
  ad::enzyme_force_circulation_gradient(
      alpha, beta, circulation_, induced_velocity_,
      bound_left_, bound_right_, control_points_, moment_arm_,
      num_vortices_, reference_area_, 0, d_circulation_enzyme);

  // Finite difference
  double d_circulation_fd[2];
  for (int i = 0; i < num_vortices_; ++i) {
    double circ_plus[2] = {circulation_[0], circulation_[1]};
    double circ_minus[2] = {circulation_[0], circulation_[1]};
    circ_plus[i] += h;
    circ_minus[i] -= h;

    double result_plus[6], result_minus[6];
    ad::compute_force_coefficients_kernel(
        alpha, beta, circ_plus, induced_velocity_,
        bound_left_, bound_right_, control_points_, moment_arm_,
        num_vortices_, reference_area_, result_plus);
    ad::compute_force_coefficients_kernel(
        alpha, beta, circ_minus, induced_velocity_,
        bound_left_, bound_right_, control_points_, moment_arm_,
        num_vortices_, reference_area_, result_minus);

    d_circulation_fd[i] = (result_plus[0] - result_minus[0]) / (2.0 * h);
  }

  for (int i = 0; i < num_vortices_; ++i) {
    EXPECT_NEAR(d_circulation_enzyme[i], d_circulation_fd[i], kTolerance)
        << "Mismatch at circulation index " << i;
  }
}

// ============================================================================
// Biot-Savart Kernel Tests
// ============================================================================

class EnzymeBiotSavartTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-6;

  void SetUp() override {
    // Set up a typical horseshoe vortex configuration
    // Control point at (0.5, 0.5, 0)
    field_point_[0] = 0.5;
    field_point_[1] = 0.5;
    field_point_[2] = 0.0;

    // Normal pointing up
    normal_[0] = 0.0;
    normal_[1] = 0.0;
    normal_[2] = 1.0;

    // Bound vortex from (0.25, 0, 0) to (0.25, 1, 0)
    bound_left_[0] = 0.25;
    bound_left_[1] = 0.0;
    bound_left_[2] = 0.0;

    bound_right_[0] = 0.25;
    bound_right_[1] = 1.0;
    bound_right_[2] = 0.0;

    prandtl_glauert_ = 1.0;  // Incompressible
    core_radius_ = 0.01;
  }

  double field_point_[3];
  double normal_[3];
  double bound_left_[3];
  double bound_right_[3];
  double prandtl_glauert_;
  double core_radius_;
};

TEST_F(EnzymeBiotSavartTest, HorseshoeVelocityIsNonZero) {
  double velocity[3];
  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, velocity);

  // Velocity should be non-zero
  double vel_mag = std::sqrt(velocity[0]*velocity[0] +
                             velocity[1]*velocity[1] +
                             velocity[2]*velocity[2]);
  EXPECT_GT(vel_mag, 0.0);
}

TEST_F(EnzymeBiotSavartTest, AICElementIsNonZero) {
  double aic = ad::compute_aic_element_kernel(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true);

  // AIC element should be non-zero for this configuration
  EXPECT_NE(aic, 0.0);
}

TEST_F(EnzymeBiotSavartTest, FieldPointGradientMatchesFD) {
  // Enzyme gradient
  double d_fp_enzyme[3];
  ad::enzyme_aic_field_point_gradient(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, d_fp_enzyme);

  // Finite difference gradient
  double d_fp_fd[3];
  ad::fd_aic_field_point_gradient(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, d_fp_fd);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(d_fp_enzyme[i], d_fp_fd[i], kTolerance)
        << "Mismatch at dimension " << i;
  }
}

TEST_F(EnzymeBiotSavartTest, BoundGradientMatchesFD) {
  // Enzyme gradient
  double d_bl_enzyme[3], d_br_enzyme[3];
  ad::enzyme_aic_bound_gradient(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, d_bl_enzyme, d_br_enzyme);

  // Finite difference gradient
  double d_bl_fd[3], d_br_fd[3];
  ad::fd_aic_bound_gradient(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, d_bl_fd, d_br_fd);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(d_bl_enzyme[i], d_bl_fd[i], kTolerance)
        << "bound_left mismatch at dimension " << i;
    EXPECT_NEAR(d_br_enzyme[i], d_br_fd[i], kTolerance)
        << "bound_right mismatch at dimension " << i;
  }
}

TEST_F(EnzymeBiotSavartTest, NormalGradientIsVelocity) {
  // The AIC element is V · n, so d(AIC)/dn = V
  double d_normal[3];
  ad::enzyme_aic_normal_gradient(
      field_point_, normal_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, d_normal);

  // Compute velocity directly
  double velocity[3];
  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, velocity);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(d_normal[i], velocity[i], kTolerance)
        << "Mismatch at dimension " << i;
  }
}

TEST_F(EnzymeBiotSavartTest, SelfInfluenceExcludesBoundVortex) {
  // When include_bound=false, only trailing legs contribute
  double vel_with_bound[3], vel_no_bound[3];

  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, true, vel_with_bound);

  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      prandtl_glauert_, core_radius_, false, vel_no_bound);

  // Should be different
  bool different = false;
  for (int i = 0; i < 3; ++i) {
    if (std::abs(vel_with_bound[i] - vel_no_bound[i]) > 1e-10) {
      different = true;
      break;
    }
  }
  EXPECT_TRUE(different);
}

TEST_F(EnzymeBiotSavartTest, CompressibilityAffectsVelocity) {
  double vel_incomp[3], vel_comp[3];

  // Incompressible
  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      1.0, core_radius_, true, vel_incomp);

  // Compressible (M=0.5, PG = sqrt(1-0.25) = 0.866)
  ad::compute_horseshoe_velocity_kernel(
      field_point_, bound_left_, bound_right_,
      0.866, core_radius_, true, vel_comp);

  // Should be different (compressibility increases velocity)
  bool different = false;
  for (int i = 0; i < 3; ++i) {
    if (std::abs(vel_incomp[i] - vel_comp[i]) > 1e-10) {
      different = true;
      break;
    }
  }
  EXPECT_TRUE(different);
}

// ============================================================================
// Mesh Kernel Tests
// ============================================================================

class EnzymeMeshKernelTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-10;

  void SetUp() override {
    // Section data for a simple tapered wing
    chord_left_ = 2.0;
    chord_right_ = 1.0;
    twist_left_ = 0.05;   // ~2.9 degrees
    twist_right_ = -0.02; // ~-1.1 degrees

    le_left_[0] = 0.0;
    le_left_[1] = 0.0;
    le_left_[2] = 0.0;

    le_right_[0] = 0.3;  // 30 deg sweep
    le_right_[1] = 5.0;  // 5m span
    le_right_[2] = 0.5;  // Some dihedral
  }

  double chord_left_, chord_right_;
  double twist_left_, twist_right_;
  double le_left_[3], le_right_[3];
};

TEST_F(EnzymeMeshKernelTest, InterpolateScalarAtEdges) {
  // At eta=0, should get left value
  EXPECT_NEAR(ad::interpolate_scalar(chord_left_, chord_right_, 0.0),
              chord_left_, kTolerance);

  // At eta=1, should get right value
  EXPECT_NEAR(ad::interpolate_scalar(chord_left_, chord_right_, 1.0),
              chord_right_, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, InterpolateScalarAtCenter) {
  // At eta=0.5, should get average
  double expected = 0.5 * (chord_left_ + chord_right_);
  EXPECT_NEAR(ad::interpolate_scalar(chord_left_, chord_right_, 0.5),
              expected, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, PanelLeadingEdgeInterpolation) {
  double le_center[3];
  ad::compute_panel_leading_edge(le_left_, le_right_, 0.5, le_center);

  for (int i = 0; i < 3; ++i) {
    double expected = 0.5 * (le_left_[i] + le_right_[i]);
    EXPECT_NEAR(le_center[i], expected, kTolerance);
  }
}

TEST_F(EnzymeMeshKernelTest, BoundVortexPositionDependsOnChord) {
  double le[3] = {0.0, 1.0, 0.0};
  double chord = 2.0;
  double x_vr = 0.25;  // Quarter chord

  double bound[3];
  ad::compute_bound_vortex_position(le, chord, x_vr, bound);

  // bound.x = le.x + x_vr * chord = 0 + 0.25 * 2 = 0.5
  EXPECT_NEAR(bound[0], 0.5, kTolerance);
  EXPECT_NEAR(bound[1], le[1], kTolerance);
  EXPECT_NEAR(bound[2], le[2], kTolerance);
}

TEST_F(EnzymeMeshKernelTest, ControlPointAt3QuarterChord) {
  double le[3] = {0.0, 1.0, 0.0};
  double chord = 2.0;
  double x_cp = 0.75;

  double cp[3];
  ad::compute_control_point_position(le, chord, x_cp, cp);

  // cp.x = le.x + x_cp * chord = 0 + 0.75 * 2 = 1.5
  EXPECT_NEAR(cp[0], 1.5, kTolerance);
  EXPECT_NEAR(cp[1], le[1], kTolerance);
  EXPECT_NEAR(cp[2], le[2], kTolerance);
}

TEST_F(EnzymeMeshKernelTest, TrefftzNormalPerpendicularToSpan) {
  double trefftz_ny, trefftz_nz;
  ad::compute_trefftz_normal(le_left_, le_right_, &trefftz_ny, &trefftz_nz);

  // Trefftz normal should be perpendicular to span direction in y-z plane
  double dy = le_right_[1] - le_left_[1];
  double dz = le_right_[2] - le_left_[2];

  // Dot product with (0, dy, dz) should be zero
  double dot = trefftz_ny * dy + trefftz_nz * dz;
  EXPECT_NEAR(dot, 0.0, kTolerance);

  // Should be unit length
  double mag = std::sqrt(trefftz_ny * trefftz_ny + trefftz_nz * trefftz_nz);
  EXPECT_NEAR(mag, 1.0, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, TwistInterpolationChordWeighted) {
  double twist = ad::interpolate_twist(chord_left_, chord_right_,
                                        twist_left_, twist_right_, 0.5);

  // Should be between left and right values
  EXPECT_GT(twist, std::min(twist_left_, twist_right_) - 0.01);
  EXPECT_LT(twist, std::max(twist_left_, twist_right_) + 0.01);
}

TEST_F(EnzymeMeshKernelTest, ChordGradientIsLinear) {
  double d_chord_left, d_chord_right;
  double eta = 0.3;

  ad::compute_chord_gradient(eta, &d_chord_left, &d_chord_right);

  // d(chord_center)/d(chord_left) = 1 - eta
  EXPECT_NEAR(d_chord_left, 1.0 - eta, kTolerance);

  // d(chord_center)/d(chord_right) = eta
  EXPECT_NEAR(d_chord_right, eta, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, LeadingEdgeGradientIsLinear) {
  double d_le_left, d_le_right;
  double eta = 0.7;

  ad::compute_leading_edge_gradient(eta, &d_le_left, &d_le_right);

  EXPECT_NEAR(d_le_left, 1.0 - eta, kTolerance);
  EXPECT_NEAR(d_le_right, eta, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, ControlPointChordGradient) {
  double x_cp = 0.75;
  double grad = ad::compute_control_point_chord_gradient(x_cp);

  // d(cp.x)/d(chord) = x_cp
  EXPECT_NEAR(grad, x_cp, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, BoundVortexChordGradient) {
  double x_vr = 0.25;
  double grad = ad::compute_bound_vortex_chord_gradient(x_vr);

  // d(bound.x)/d(chord) = x_vr
  EXPECT_NEAR(grad, x_vr, kTolerance);
}

TEST_F(EnzymeMeshKernelTest, PanelNormalIsUnitVector) {
  double bound_left[3] = {0.5, 0.0, 0.0};
  double bound_right[3] = {0.5, 5.0, 0.5};
  double twist = 0.0;
  double trefftz_ny, trefftz_nz;
  ad::compute_trefftz_normal(le_left_, le_right_, &trefftz_ny, &trefftz_nz);

  double normal[3];
  ad::compute_panel_normal(twist, trefftz_ny, trefftz_nz,
                           bound_left, bound_right, 0.0, normal);

  double mag = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  EXPECT_NEAR(mag, 1.0, kTolerance);
}

// ============================================================================
// Trefftz Plane Kernel Tests
// ============================================================================

class EnzymeTrefftzTest : public ::testing::Test {
 protected:
  static constexpr double kTolerance = 1e-6;

  void SetUp() override {
    // Simple 3-strip wing (symmetric, planar)
    num_strips_ = 3;

    // Strip circulations (typical elliptic distribution)
    strip_gamma_[0] = 1.0;
    strip_gamma_[1] = 1.2;
    strip_gamma_[2] = 1.0;

    // Trailing edge positions (left endpoints)
    // Strip 0: y = -5 to -1.67
    strip_te_left_[0] = 1.0; strip_te_left_[1] = -5.0; strip_te_left_[2] = 0.0;
    strip_te_left_[3] = 1.0; strip_te_left_[4] = -1.67; strip_te_left_[5] = 0.0;
    strip_te_left_[6] = 1.0; strip_te_left_[7] = 1.67; strip_te_left_[8] = 0.0;

    // Trailing edge positions (right endpoints)
    strip_te_right_[0] = 1.0; strip_te_right_[1] = -1.67; strip_te_right_[2] = 0.0;
    strip_te_right_[3] = 1.0; strip_te_right_[4] = 1.67; strip_te_right_[5] = 0.0;
    strip_te_right_[6] = 1.0; strip_te_right_[7] = 5.0; strip_te_right_[8] = 0.0;

    // All strips on same surface (index 0)
    strip_surface_index_[0] = 0;
    strip_surface_index_[1] = 0;
    strip_surface_index_[2] = 0;

    // Chord lengths
    strip_chord_[0] = 1.0;
    strip_chord_[1] = 1.0;
    strip_chord_[2] = 1.0;

    reference_area_ = 10.0;
    reference_span_ = 10.0;
    core_chord_fraction_ = 0.25;
    core_width_fraction_ = 0.25;
  }

  int num_strips_;
  double strip_gamma_[3];
  double strip_te_left_[9];
  double strip_te_right_[9];
  int strip_surface_index_[3];
  double strip_chord_[3];
  double reference_area_;
  double reference_span_;
  double core_chord_fraction_;
  double core_width_fraction_;
};

TEST_F(EnzymeTrefftzTest, TrefftzDragIsNonNegative) {
  double CL_ff, CY_ff, CDi_ff;
  ad::compute_trefftz_drag_kernel(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      &CL_ff, &CY_ff, &CDi_ff);

  // Induced drag should be non-negative for reasonable configurations
  EXPECT_GE(CDi_ff, 0.0);
}

TEST_F(EnzymeTrefftzTest, TrefftzLiftIsPositive) {
  double CL_ff, CY_ff, CDi_ff;
  ad::compute_trefftz_drag_kernel(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      &CL_ff, &CY_ff, &CDi_ff);

  // Positive circulation should give positive lift
  EXPECT_GT(CL_ff, 0.0);
}

TEST_F(EnzymeTrefftzTest, TrefftzSideForceIsSmall) {
  double CL_ff, CY_ff, CDi_ff;
  ad::compute_trefftz_drag_kernel(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      &CL_ff, &CY_ff, &CDi_ff);

  // For planar symmetric wing, side force should be zero
  EXPECT_NEAR(CY_ff, 0.0, kTolerance);
}

TEST_F(EnzymeTrefftzTest, SpanEfficiencyIsReasonable) {
  double CL_ff, CY_ff, CDi_ff;
  ad::compute_trefftz_drag_kernel(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      &CL_ff, &CY_ff, &CDi_ff);

  double e = ad::compute_span_efficiency(CL_ff, CY_ff, CDi_ff,
                                          reference_area_, reference_span_);

  // Span efficiency should be between 0 and slightly above 1
  // (can exceed 1 for non-planar configurations, but should be reasonable)
  EXPECT_GT(e, 0.0);
  EXPECT_LT(e, 2.0);  // Reasonable upper bound
}

TEST_F(EnzymeTrefftzTest, CirculationGradientMatchesFD) {
  // Test gradient of CDi_ff (output_index=2) w.r.t. circulation
  double d_gamma_enzyme[3];
  ad::enzyme_trefftz_circulation_gradient(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      2,  // CDi_ff
      d_gamma_enzyme);

  double d_gamma_fd[3];
  ad::fd_trefftz_circulation_gradient(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      2,  // CDi_ff
      d_gamma_fd);

  for (int i = 0; i < num_strips_; ++i) {
    EXPECT_NEAR(d_gamma_enzyme[i], d_gamma_fd[i], kTolerance)
        << "Mismatch at strip " << i;
  }
}

TEST_F(EnzymeTrefftzTest, CLGradientMatchesFD) {
  // Test gradient of CL_ff (output_index=0) w.r.t. circulation
  double d_gamma_enzyme[3];
  ad::enzyme_trefftz_circulation_gradient(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      0,  // CL_ff
      d_gamma_enzyme);

  double d_gamma_fd[3];
  ad::fd_trefftz_circulation_gradient(
      strip_gamma_, strip_te_left_, strip_te_right_, strip_chord_,
      strip_surface_index_, num_strips_, reference_area_, reference_span_,
      core_chord_fraction_, core_width_fraction_,
      0,  // CL_ff
      d_gamma_fd);

  for (int i = 0; i < num_strips_; ++i) {
    EXPECT_NEAR(d_gamma_enzyme[i], d_gamma_fd[i], kTolerance)
        << "Mismatch at strip " << i;
  }
}

#else  // !VLM_USE_ENZYME

// When Enzyme is not available, skip all tests gracefully

TEST_F(EnzymeBasicTest, SkippedWithoutEnzyme) {
  GTEST_SKIP() << "Enzyme not available - enable with -DVLM_ENABLE_ENZYME=ON";
}

#endif  // VLM_USE_ENZYME

}  // namespace
}  // namespace vlm
