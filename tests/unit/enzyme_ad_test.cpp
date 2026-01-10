// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <cmath>

#include <gtest/gtest.h>

#include "vlm/ad/enzyme_interface.h"

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

#else  // !VLM_USE_ENZYME

// When Enzyme is not available, skip all tests gracefully

TEST_F(EnzymeBasicTest, SkippedWithoutEnzyme) {
  GTEST_SKIP() << "Enzyme not available - enable with -DVLM_ENABLE_ENZYME=ON";
}

#endif  // VLM_USE_ENZYME

}  // namespace
}  // namespace vlm
