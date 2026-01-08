// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <cmath>
#include <functional>

#include <gtest/gtest.h>

#include "vlm/ad/autodiff_tape.h"
#include "vlm/ad/dual_number.h"

namespace vlm {
namespace ad {
namespace {

constexpr double kTestTolerance = 1.0e-8;
constexpr double kFiniteDiffStep = 1.0e-6;
constexpr double kFiniteDiffTolerance = 1.0e-5;

// Helper: compute finite difference derivative
double FiniteDifference(std::function<double(double)> f, double x,
                        double h = kFiniteDiffStep) {
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

// Helper: compute finite difference gradient
std::vector<double> FiniteDifferenceGradient(
    std::function<double(const std::vector<double>&)> f,
    const std::vector<double>& x, double h = kFiniteDiffStep) {
  std::vector<double> grad(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    std::vector<double> x_plus = x;
    std::vector<double> x_minus = x;
    x_plus[i] += h;
    x_minus[i] -= h;
    grad[i] = (f(x_plus) - f(x_minus)) / (2.0 * h);
  }
  return grad;
}

// ============================================================================
// Dual Number Tests - Forward Mode AD
// ============================================================================

class DualNumberTest : public ::testing::Test {};

TEST_F(DualNumberTest, DefaultConstruction) {
  Dual1d x;
  EXPECT_DOUBLE_EQ(x.Value(), 0.0);
  EXPECT_DOUBLE_EQ(x.Derivative(), 0.0);
}

TEST_F(DualNumberTest, ValueOnlyConstruction) {
  Dual1d x(3.0);
  EXPECT_DOUBLE_EQ(x.Value(), 3.0);
  EXPECT_DOUBLE_EQ(x.Derivative(), 0.0);
}

TEST_F(DualNumberTest, FullConstruction) {
  Dual1d x(3.0, 1.0);
  EXPECT_DOUBLE_EQ(x.Value(), 3.0);
  EXPECT_DOUBLE_EQ(x.Derivative(), 1.0);
}

TEST_F(DualNumberTest, VariableFactory) {
  auto x = Dual1d::Variable(2.0);
  EXPECT_DOUBLE_EQ(x.Value(), 2.0);
  EXPECT_DOUBLE_EQ(x.Derivative(), 1.0);
}

TEST_F(DualNumberTest, ConstantFactory) {
  auto c = Dual1d::Constant(5.0);
  EXPECT_DOUBLE_EQ(c.Value(), 5.0);
  EXPECT_DOUBLE_EQ(c.Derivative(), 0.0);
}

TEST_F(DualNumberTest, Addition) {
  auto x = Dual1d::Variable(2.0);
  Dual1d c(3.0);
  auto y = x + c;

  EXPECT_DOUBLE_EQ(y.Value(), 5.0);
  EXPECT_DOUBLE_EQ(y.Derivative(), 1.0);
}

TEST_F(DualNumberTest, Subtraction) {
  auto x = Dual1d::Variable(5.0);
  Dual1d c(3.0);
  auto y = x - c;

  EXPECT_DOUBLE_EQ(y.Value(), 2.0);
  EXPECT_DOUBLE_EQ(y.Derivative(), 1.0);
}

TEST_F(DualNumberTest, Multiplication) {
  auto x = Dual1d::Variable(3.0);
  Dual1d c(4.0);
  auto y = x * c;  // y = 4x

  EXPECT_DOUBLE_EQ(y.Value(), 12.0);
  EXPECT_DOUBLE_EQ(y.Derivative(), 4.0);
}

TEST_F(DualNumberTest, Division) {
  auto x = Dual1d::Variable(6.0);
  Dual1d c(2.0);
  auto y = x / c;  // y = x/2

  EXPECT_DOUBLE_EQ(y.Value(), 3.0);
  EXPECT_DOUBLE_EQ(y.Derivative(), 0.5);
}

TEST_F(DualNumberTest, Negation) {
  auto x = Dual1d::Variable(3.0);
  auto y = -x;

  EXPECT_DOUBLE_EQ(y.Value(), -3.0);
  EXPECT_DOUBLE_EQ(y.Derivative(), -1.0);
}

TEST_F(DualNumberTest, Polynomial) {
  // f(x) = x³ + 2x² - 5x + 3
  // f'(x) = 3x² + 4x - 5
  auto x = Dual1d::Variable(2.0);
  auto f = x * x * x + 2.0 * x * x - 5.0 * x + 3.0;

  double expected_value = 8.0 + 8.0 - 10.0 + 3.0;  // 9
  double expected_deriv = 12.0 + 8.0 - 5.0;        // 15

  EXPECT_NEAR(f.Value(), expected_value, kTestTolerance);
  EXPECT_NEAR(f.Derivative(), expected_deriv, kTestTolerance);
}

TEST_F(DualNumberTest, SinCos) {
  double x0 = 0.5;
  auto x = Dual1d::Variable(x0);

  auto s = sin(x);
  auto c = cos(x);

  EXPECT_NEAR(s.Value(), std::sin(x0), kTestTolerance);
  EXPECT_NEAR(s.Derivative(), std::cos(x0), kTestTolerance);  // d/dx sin(x) = cos(x)

  EXPECT_NEAR(c.Value(), std::cos(x0), kTestTolerance);
  EXPECT_NEAR(c.Derivative(), -std::sin(x0), kTestTolerance);  // d/dx cos(x) = -sin(x)
}

TEST_F(DualNumberTest, ExpLog) {
  double x0 = 2.0;
  auto x = Dual1d::Variable(x0);

  auto e = exp(x);
  auto l = log(x);

  EXPECT_NEAR(e.Value(), std::exp(x0), kTestTolerance);
  EXPECT_NEAR(e.Derivative(), std::exp(x0), kTestTolerance);  // d/dx exp(x) = exp(x)

  EXPECT_NEAR(l.Value(), std::log(x0), kTestTolerance);
  EXPECT_NEAR(l.Derivative(), 1.0 / x0, kTestTolerance);  // d/dx log(x) = 1/x
}

TEST_F(DualNumberTest, Sqrt) {
  double x0 = 4.0;
  auto x = Dual1d::Variable(x0);

  auto s = sqrt(x);

  EXPECT_NEAR(s.Value(), 2.0, kTestTolerance);
  EXPECT_NEAR(s.Derivative(), 0.25, kTestTolerance);  // d/dx sqrt(x) = 0.5/sqrt(x)
}

TEST_F(DualNumberTest, PowerConstantExponent) {
  double x0 = 2.0;
  double n = 3.0;
  auto x = Dual1d::Variable(x0);

  auto p = pow(x, n);

  EXPECT_NEAR(p.Value(), 8.0, kTestTolerance);
  EXPECT_NEAR(p.Derivative(), 12.0, kTestTolerance);  // d/dx x³ = 3x²
}

TEST_F(DualNumberTest, PowerConstantBase) {
  double a = 2.0;
  double x0 = 3.0;
  auto x = Dual1d::Variable(x0);

  auto p = pow(a, x);  // 2^x

  double expected_val = std::pow(a, x0);
  double expected_deriv = expected_val * std::log(a);

  EXPECT_NEAR(p.Value(), expected_val, kTestTolerance);
  EXPECT_NEAR(p.Derivative(), expected_deriv, kTestTolerance);
}

TEST_F(DualNumberTest, ChainRule) {
  // f(x) = sin(x²)
  // f'(x) = 2x * cos(x²)
  double x0 = 1.5;
  auto x = Dual1d::Variable(x0);

  auto f = sin(x * x);

  double expected_val = std::sin(x0 * x0);
  double expected_deriv = 2.0 * x0 * std::cos(x0 * x0);

  EXPECT_NEAR(f.Value(), expected_val, kTestTolerance);
  EXPECT_NEAR(f.Derivative(), expected_deriv, kTestTolerance);
}

TEST_F(DualNumberTest, ComplexExpression) {
  // f(x) = exp(-x²/2) / sqrt(2*pi)
  // f'(x) = -x * exp(-x²/2) / sqrt(2*pi)
  double x0 = 1.0;
  auto x = Dual1d::Variable(x0);

  constexpr double two_pi = 2.0 * 3.14159265358979323846;
  auto f = exp(-x * x / 2.0) / sqrt(Dual1d(two_pi));

  double expected_val = std::exp(-x0 * x0 / 2.0) / std::sqrt(two_pi);
  double expected_deriv = -x0 * std::exp(-x0 * x0 / 2.0) / std::sqrt(two_pi);

  EXPECT_NEAR(f.Value(), expected_val, kTestTolerance);
  EXPECT_NEAR(f.Derivative(), expected_deriv, kTestTolerance);
}

TEST_F(DualNumberTest, MultiVariableGradient) {
  // f(a, b) = a² + a*b + b²
  // ∂f/∂a = 2a + b, ∂f/∂b = a + 2b
  double a0 = 2.0, b0 = 3.0;

  auto a = Dual<double, 2>::Variable(a0, 0);  // ∂/∂a
  auto b = Dual<double, 2>::Variable(b0, 1);  // ∂/∂b

  auto f = a * a + a * b + b * b;

  double expected_val = a0 * a0 + a0 * b0 + b0 * b0;  // 4 + 6 + 9 = 19
  double expected_da = 2.0 * a0 + b0;                  // 4 + 3 = 7
  double expected_db = a0 + 2.0 * b0;                  // 2 + 6 = 8

  EXPECT_NEAR(f.Value(), expected_val, kTestTolerance);
  EXPECT_NEAR(f.Derivative(0), expected_da, kTestTolerance);
  EXPECT_NEAR(f.Derivative(1), expected_db, kTestTolerance);
}

TEST_F(DualNumberTest, ComparisonOperators) {
  Dual1d a(3.0, 1.0);
  Dual1d b(5.0, 2.0);
  Dual1d c(3.0, 3.0);

  EXPECT_TRUE(a < b);
  EXPECT_TRUE(a <= c);
  EXPECT_TRUE(b > a);
  EXPECT_TRUE(a == c);  // Compares values only
  EXPECT_TRUE(a != b);
}

TEST_F(DualNumberTest, VerifyAgainstFiniteDifference) {
  // Test several functions against finite differences
  auto test_func = [](std::function<Dual1d(Dual1d)> f_dual,
                      std::function<double(double)> f_double, double x0) {
    auto x = Dual1d::Variable(x0);
    auto y = f_dual(x);
    double fd_deriv = FiniteDifference(f_double, x0);
    EXPECT_NEAR(y.Derivative(), fd_deriv, kFiniteDiffTolerance)
        << "Mismatch at x = " << x0;
  };

  // f(x) = x³ - 2x + 1
  test_func([](Dual1d x) { return x * x * x - 2.0 * x + 1.0; },
            [](double x) { return x * x * x - 2.0 * x + 1.0; }, 1.5);

  // f(x) = sin(x) * cos(x)
  test_func([](Dual1d x) { return sin(x) * cos(x); },
            [](double x) { return std::sin(x) * std::cos(x); }, 0.7);

  // f(x) = exp(-x²)
  test_func([](Dual1d x) { return exp(-x * x); },
            [](double x) { return std::exp(-x * x); }, 0.5);

  // f(x) = sqrt(1 + x²)
  test_func([](Dual1d x) { return sqrt(1.0 + x * x); },
            [](double x) { return std::sqrt(1.0 + x * x); }, 2.0);

  // f(x) = log(1 + exp(x))  (softplus)
  test_func([](Dual1d x) { return log(1.0 + exp(x)); },
            [](double x) { return std::log(1.0 + std::exp(x)); }, 1.0);
}

// ============================================================================
// Tape-Based Reverse Mode AD Tests
// ============================================================================

class TapeADTest : public ::testing::Test {};

TEST_F(TapeADTest, SimpleAddition) {
  Tape tape;
  TapeScope scope(tape);

  auto x = tape.CreateVariable(2.0);
  auto y = tape.CreateVariable(3.0);
  auto z = x + y;

  auto grad = tape.ComputeGradient(z);

  EXPECT_DOUBLE_EQ(z.Value(), 5.0);
  EXPECT_DOUBLE_EQ(grad[0], 1.0);  // dz/dx
  EXPECT_DOUBLE_EQ(grad[1], 1.0);  // dz/dy
}

TEST_F(TapeADTest, SimpleMultiplication) {
  Tape tape;
  TapeScope scope(tape);

  auto x = tape.CreateVariable(2.0);
  auto y = tape.CreateVariable(3.0);
  auto z = x * y;

  auto grad = tape.ComputeGradient(z);

  EXPECT_DOUBLE_EQ(z.Value(), 6.0);
  EXPECT_DOUBLE_EQ(grad[0], 3.0);  // dz/dx = y
  EXPECT_DOUBLE_EQ(grad[1], 2.0);  // dz/dy = x
}

TEST_F(TapeADTest, SimpleDivision) {
  Tape tape;
  TapeScope scope(tape);

  auto x = tape.CreateVariable(6.0);
  auto y = tape.CreateVariable(2.0);
  auto z = x / y;  // z = x/y = 3

  auto grad = tape.ComputeGradient(z);

  EXPECT_DOUBLE_EQ(z.Value(), 3.0);
  EXPECT_DOUBLE_EQ(grad[0], 0.5);   // dz/dx = 1/y
  EXPECT_DOUBLE_EQ(grad[1], -1.5);  // dz/dy = -x/y²
}

TEST_F(TapeADTest, Polynomial) {
  // f(x, y) = x² + 2xy + y²
  // ∂f/∂x = 2x + 2y, ∂f/∂y = 2x + 2y
  Tape tape;
  TapeScope scope(tape);

  auto x = tape.CreateVariable(3.0);
  auto y = tape.CreateVariable(4.0);
  auto f = x * x + 2.0 * x * y + y * y;

  auto grad = tape.ComputeGradient(f);

  EXPECT_DOUBLE_EQ(f.Value(), 49.0);       // (3+4)² = 49
  EXPECT_DOUBLE_EQ(grad[0], 14.0);         // 2*3 + 2*4 = 14
  EXPECT_DOUBLE_EQ(grad[1], 14.0);         // 2*3 + 2*4 = 14
}

TEST_F(TapeADTest, SinCos) {
  Tape tape;
  TapeScope scope(tape);

  double x0 = 0.5;
  auto x = tape.CreateVariable(x0);
  auto s = sin(x);
  auto c = cos(x);

  auto grad_sin = tape.ComputeGradient(s);
  auto grad_cos = tape.ComputeGradient(c);

  EXPECT_NEAR(s.Value(), std::sin(x0), kTestTolerance);
  EXPECT_NEAR(grad_sin[0], std::cos(x0), kTestTolerance);

  EXPECT_NEAR(c.Value(), std::cos(x0), kTestTolerance);
  EXPECT_NEAR(grad_cos[0], -std::sin(x0), kTestTolerance);
}

TEST_F(TapeADTest, ExpLog) {
  Tape tape;
  TapeScope scope(tape);

  double x0 = 2.0;
  auto x = tape.CreateVariable(x0);
  auto e = exp(x);
  auto l = log(x);

  auto grad_exp = tape.ComputeGradient(e);
  auto grad_log = tape.ComputeGradient(l);

  EXPECT_NEAR(e.Value(), std::exp(x0), kTestTolerance);
  EXPECT_NEAR(grad_exp[0], std::exp(x0), kTestTolerance);

  EXPECT_NEAR(l.Value(), std::log(x0), kTestTolerance);
  EXPECT_NEAR(grad_log[0], 1.0 / x0, kTestTolerance);
}

TEST_F(TapeADTest, ChainRule) {
  // f(x) = sin(x²)
  Tape tape;
  TapeScope scope(tape);

  double x0 = 1.5;
  auto x = tape.CreateVariable(x0);
  auto f = sin(x * x);

  auto grad = tape.ComputeGradient(f);

  double expected_deriv = 2.0 * x0 * std::cos(x0 * x0);
  EXPECT_NEAR(f.Value(), std::sin(x0 * x0), kTestTolerance);
  EXPECT_NEAR(grad[0], expected_deriv, kTestTolerance);
}

TEST_F(TapeADTest, RosenbrockFunction) {
  // f(x, y) = (1-x)² + 100*(y-x²)²
  // ∂f/∂x = -2(1-x) - 400*x*(y-x²)
  // ∂f/∂y = 200*(y-x²)
  Tape tape;
  TapeScope scope(tape);

  double x0 = 1.5, y0 = 2.25;
  auto x = tape.CreateVariable(x0);
  auto y = tape.CreateVariable(y0);

  auto term1 = (1.0 - x) * (1.0 - x);
  auto term2 = 100.0 * (y - x * x) * (y - x * x);
  auto f = term1 + term2;

  auto grad = tape.ComputeGradient(f);

  double diff = y0 - x0 * x0;
  double expected_dx = -2.0 * (1.0 - x0) - 400.0 * x0 * diff;
  double expected_dy = 200.0 * diff;

  EXPECT_NEAR(grad[0], expected_dx, kTestTolerance);
  EXPECT_NEAR(grad[1], expected_dy, kTestTolerance);
}

TEST_F(TapeADTest, ComputeGradientHelper) {
  // Use the helper function
  auto f = [](const std::vector<TapeVariable>& x) {
    return x[0] * x[0] + x[1] * x[1] + x[0] * x[1];
  };

  std::vector<double> point = {2.0, 3.0};
  auto grad = ComputeGradient(f, point);

  // f(x,y) = x² + y² + xy
  // ∂f/∂x = 2x + y = 7
  // ∂f/∂y = 2y + x = 8
  EXPECT_NEAR(grad[0], 7.0, kTestTolerance);
  EXPECT_NEAR(grad[1], 8.0, kTestTolerance);
}

TEST_F(TapeADTest, VerifyAgainstFiniteDifference) {
  // Test gradient computation against finite differences
  auto f_tape = [](const std::vector<TapeVariable>& x) {
    return exp(-x[0] * x[0] - x[1] * x[1]) * sin(x[0] + x[1]);
  };

  auto f_double = [](const std::vector<double>& x) {
    return std::exp(-x[0] * x[0] - x[1] * x[1]) * std::sin(x[0] + x[1]);
  };

  std::vector<double> point = {0.5, 0.7};

  auto ad_grad = ComputeGradient(f_tape, point);
  auto fd_grad = FiniteDifferenceGradient(f_double, point);

  EXPECT_NEAR(ad_grad[0], fd_grad[0], kFiniteDiffTolerance);
  EXPECT_NEAR(ad_grad[1], fd_grad[1], kFiniteDiffTolerance);
}

TEST_F(TapeADTest, ManyVariables) {
  // Test with more variables
  int n = 10;
  Tape tape;
  TapeScope scope(tape);

  std::vector<TapeVariable> x;
  for (int i = 0; i < n; ++i) {
    x.push_back(tape.CreateVariable(static_cast<double>(i + 1)));
  }

  // f = sum of x_i²
  TapeVariable f = x[0] * x[0];
  for (int i = 1; i < n; ++i) {
    f = f + x[i] * x[i];
  }

  auto grad = tape.ComputeGradient(f);

  // ∂f/∂x_i = 2*x_i
  for (int i = 0; i < n; ++i) {
    double expected = 2.0 * (i + 1);
    EXPECT_NEAR(grad[i], expected, kTestTolerance);
  }
}

TEST_F(TapeADTest, TrigFunctions) {
  Tape tape;
  TapeScope scope(tape);

  double x0 = 0.3;
  auto x = tape.CreateVariable(x0);

  auto t = tan(x);
  auto as = asin(x);
  auto ac = acos(x);
  auto at = atan(x);

  // tan
  EXPECT_NEAR(t.Value(), std::tan(x0), kTestTolerance);
  auto grad_tan = tape.ComputeGradient(t);
  double c = std::cos(x0);
  EXPECT_NEAR(grad_tan[0], 1.0 / (c * c), kTestTolerance);

  // asin
  EXPECT_NEAR(as.Value(), std::asin(x0), kTestTolerance);
  auto grad_asin = tape.ComputeGradient(as);
  EXPECT_NEAR(grad_asin[0], 1.0 / std::sqrt(1.0 - x0 * x0), kTestTolerance);

  // acos
  EXPECT_NEAR(ac.Value(), std::acos(x0), kTestTolerance);
  auto grad_acos = tape.ComputeGradient(ac);
  EXPECT_NEAR(grad_acos[0], -1.0 / std::sqrt(1.0 - x0 * x0), kTestTolerance);

  // atan
  EXPECT_NEAR(at.Value(), std::atan(x0), kTestTolerance);
  auto grad_atan = tape.ComputeGradient(at);
  EXPECT_NEAR(grad_atan[0], 1.0 / (1.0 + x0 * x0), kTestTolerance);
}

TEST_F(TapeADTest, HyperbolicFunctions) {
  Tape tape;
  TapeScope scope(tape);

  double x0 = 0.5;
  auto x = tape.CreateVariable(x0);

  auto sh = sinh(x);
  auto ch = cosh(x);
  auto th = tanh(x);

  EXPECT_NEAR(sh.Value(), std::sinh(x0), kTestTolerance);
  EXPECT_NEAR(ch.Value(), std::cosh(x0), kTestTolerance);
  EXPECT_NEAR(th.Value(), std::tanh(x0), kTestTolerance);

  auto grad_sinh = tape.ComputeGradient(sh);
  auto grad_cosh = tape.ComputeGradient(ch);
  auto grad_tanh = tape.ComputeGradient(th);

  EXPECT_NEAR(grad_sinh[0], std::cosh(x0), kTestTolerance);
  EXPECT_NEAR(grad_cosh[0], std::sinh(x0), kTestTolerance);
  double th_val = std::tanh(x0);
  EXPECT_NEAR(grad_tanh[0], 1.0 - th_val * th_val, kTestTolerance);
}

TEST_F(TapeADTest, ScalarOperations) {
  Tape tape;
  TapeScope scope(tape);

  auto x = tape.CreateVariable(3.0);

  auto a = 2.0 + x;  // = 5
  auto b = 2.0 - x;  // = -1
  auto c = 2.0 * x;  // = 6
  auto d = 6.0 / x;  // = 2

  EXPECT_DOUBLE_EQ(a.Value(), 5.0);
  EXPECT_DOUBLE_EQ(b.Value(), -1.0);
  EXPECT_DOUBLE_EQ(c.Value(), 6.0);
  EXPECT_DOUBLE_EQ(d.Value(), 2.0);

  auto grad_a = tape.ComputeGradient(a);
  auto grad_b = tape.ComputeGradient(b);
  auto grad_c = tape.ComputeGradient(c);
  auto grad_d = tape.ComputeGradient(d);

  EXPECT_DOUBLE_EQ(grad_a[0], 1.0);
  EXPECT_DOUBLE_EQ(grad_b[0], -1.0);
  EXPECT_DOUBLE_EQ(grad_c[0], 2.0);
  EXPECT_NEAR(grad_d[0], -6.0 / 9.0, kTestTolerance);
}

// ============================================================================
// Cross-Validation: Forward vs Reverse Mode
// ============================================================================

TEST(CrossValidationTest, ForwardVsReverseSimple) {
  // f(x) = x³ + sin(x)
  double x0 = 2.0;

  // Forward mode
  auto xf = Dual1d::Variable(x0);
  auto ff = xf * xf * xf + sin(xf);
  double forward_deriv = ff.Derivative();

  // Reverse mode
  auto f_tape = [](const std::vector<TapeVariable>& x) {
    return x[0] * x[0] * x[0] + sin(x[0]);
  };
  auto reverse_grad = ComputeGradient(f_tape, {x0});

  EXPECT_NEAR(forward_deriv, reverse_grad[0], kTestTolerance);

  // Also verify against finite difference
  auto f_double = [](double x) {
    return x * x * x + std::sin(x);
  };
  double fd_deriv = FiniteDifference(f_double, x0);

  EXPECT_NEAR(forward_deriv, fd_deriv, kFiniteDiffTolerance);
}

TEST(CrossValidationTest, ForwardVsReverseMultiVar) {
  // f(x, y) = exp(x*y) + x²*y
  double x0 = 1.0, y0 = 2.0;

  // Forward mode for ∂f/∂x
  {
    Dual<double, 2> x(x0, {1.0, 0.0});
    Dual<double, 2> y(y0, {0.0, 1.0});
    auto f = exp(x * y) + x * x * y;

    // Reverse mode
    auto f_tape = [](const std::vector<TapeVariable>& v) {
      return exp(v[0] * v[1]) + v[0] * v[0] * v[1];
    };
    auto reverse_grad = ComputeGradient(f_tape, {x0, y0});

    EXPECT_NEAR(f.Derivative(0), reverse_grad[0], kTestTolerance);
    EXPECT_NEAR(f.Derivative(1), reverse_grad[1], kTestTolerance);
  }
}

}  // namespace
}  // namespace ad
}  // namespace vlm
