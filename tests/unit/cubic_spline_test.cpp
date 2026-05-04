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

#include "pegasus/utils/cubic_spline.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

// Test tolerance for floating point comparisons
constexpr double kTestTolerance = 1.0e-10;

class CubicSplineTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

TEST_F(CubicSplineTest, LinearDataInterpolation) {
    // A cubic spline through linear data should reproduce the line exactly
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<Scalar> x = {0.0, 2.0, 4.0, 6.0, 8.0};  // x = 2*s

    CubicSpline spline;
    spline.Build(s, x);

    // Test interpolation at various points
    EXPECT_NEAR(spline.Evaluate(0.5), 1.0, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(1.5), 3.0, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(2.5), 5.0, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(3.5), 7.0, kTestTolerance);

    // Derivative should be constant = 2
    EXPECT_NEAR(spline.EvaluateDerivative(0.5), 2.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(2.0), 2.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(3.5), 2.0, kTestTolerance);

    // Second derivative should be zero for linear data
    EXPECT_NEAR(spline.EvaluateSecondDerivative(1.0), 0.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateSecondDerivative(2.5), 0.0, kTestTolerance);
}

TEST_F(CubicSplineTest, QuadraticDataWithMatchingBoundaryConditions) {
    // A cubic spline can exactly reproduce a quadratic if we specify
    // the correct boundary slopes (derivative of x = s^2 is 2*s)
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<Scalar> x(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        x[i] = s[i] * s[i];  // x = s^2
    }

    CubicSpline spline;
    // Specify the exact slopes at endpoints: dx/ds = 2*s
    // At s=0: slope = 0, at s=4: slope = 8
    spline.Build(s, x, SplineBoundaryCondition::kSpecifiedSlope,
                 SplineBoundaryCondition::kSpecifiedSlope, 0.0, 8.0);

    // Test interpolation - should be exact with correct BCs
    EXPECT_NEAR(spline.Evaluate(0.5), 0.25, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(1.5), 2.25, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(2.5), 6.25, kTestTolerance);

    // Derivative should be 2*s
    EXPECT_NEAR(spline.EvaluateDerivative(0.5), 1.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(1.5), 3.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(2.5), 5.0, kTestTolerance);

    // Second derivative should be constant = 2
    EXPECT_NEAR(spline.EvaluateSecondDerivative(1.0), 2.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateSecondDerivative(2.5), 2.0, kTestTolerance);
}

TEST_F(CubicSplineTest, SinusoidalData) {
    // Test with sinusoidal data - should be smooth and reasonably accurate
    const int num_points = 21;
    std::vector<Scalar> s(num_points);
    std::vector<Scalar> x(num_points);

    for (int i = 0; i < num_points; ++i) {
        s[i] = 2.0 * kPi * static_cast<Scalar>(i) / (num_points - 1);
        x[i] = std::sin(s[i]);
    }

    CubicSpline spline;
    spline.Build(s, x);

    // Test at midpoints (more challenging for interpolation)
    for (int i = 0; i < num_points - 1; ++i) {
        Scalar s_mid = 0.5 * (s[i] + s[i + 1]);
        Scalar expected = std::sin(s_mid);
        Scalar actual = spline.Evaluate(s_mid);
        // Allow slightly larger tolerance for transcendental functions
        EXPECT_NEAR(actual, expected, 1.0e-4);
    }

    // Test derivative at a known point
    Scalar s_test = kPi / 4.0;
    Scalar expected_deriv = std::cos(s_test);
    Scalar actual_deriv = spline.EvaluateDerivative(s_test);
    EXPECT_NEAR(actual_deriv, expected_deriv, 1.0e-3);
}

TEST_F(CubicSplineTest, EvaluateAllConsistency) {
    // Verify that EvaluateAll returns the same results as individual calls
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0};
    std::vector<Scalar> x = {1.0, 2.0, 0.5, 3.0};

    CubicSpline spline;
    spline.Build(s, x);

    Scalar s_eval = 1.7;
    Scalar val, deriv1, deriv2;
    spline.EvaluateAll(s_eval, val, deriv1, deriv2);

    EXPECT_NEAR(val, spline.Evaluate(s_eval), kTestTolerance);
    EXPECT_NEAR(deriv1, spline.EvaluateDerivative(s_eval), kTestTolerance);
    EXPECT_NEAR(deriv2, spline.EvaluateSecondDerivative(s_eval),
                kTestTolerance);
}

TEST_F(CubicSplineTest, EndpointInterpolation) {
    // Spline should pass exactly through data points
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<Scalar> x = {1.5, 2.3, 0.7, 4.1, 2.9};

    CubicSpline spline;
    spline.Build(s, x);

    for (size_t i = 0; i < s.size(); ++i) {
        EXPECT_NEAR(spline.Evaluate(s[i]), x[i], kTestTolerance);
    }
}

TEST_F(CubicSplineTest, SpecifiedSlopeBoundaryCondition) {
    // Test with specified slopes at endpoints
    std::vector<Scalar> s = {0.0, 1.0, 2.0};
    std::vector<Scalar> x = {0.0, 1.0, 0.0};  // Parabola-like

    CubicSpline spline;
    spline.Build(s, x, SplineBoundaryCondition::kSpecifiedSlope,
                 SplineBoundaryCondition::kSpecifiedSlope, 2.0, -2.0);

    // Verify slopes at endpoints
    EXPECT_NEAR(spline.EvaluateDerivative(0.0), 2.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(2.0), -2.0, kTestTolerance);

    // Should still pass through data points
    EXPECT_NEAR(spline.Evaluate(0.0), 0.0, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(1.0), 1.0, kTestTolerance);
    EXPECT_NEAR(spline.Evaluate(2.0), 0.0, kTestTolerance);
}

TEST_F(CubicSplineTest, TwoPointSpline) {
    // Minimal case: 2 points
    std::vector<Scalar> s = {0.0, 1.0};
    std::vector<Scalar> x = {0.0, 2.0};

    CubicSpline spline;
    spline.Build(s, x);

    // Should interpolate linearly for natural BC with 2 points
    EXPECT_NEAR(spline.Evaluate(0.5), 1.0, kTestTolerance);
    EXPECT_NEAR(spline.EvaluateDerivative(0.5), 2.0, kTestTolerance);
}

TEST_F(CubicSplineTest, NonUniformSpacing) {
    // Test with non-uniformly spaced data
    std::vector<Scalar> s = {0.0, 0.1, 0.5, 2.0, 4.0};
    std::vector<Scalar> x = {1.0, 1.2, 2.0, 3.0, 2.5};

    CubicSpline spline;
    spline.Build(s, x);

    // Verify interpolation at data points
    for (size_t i = 0; i < s.size(); ++i) {
        EXPECT_NEAR(spline.Evaluate(s[i]), x[i], kTestTolerance);
    }

    // Verify smoothness - second derivative should be continuous
    // Check at a point near an interior knot
    Scalar s_left = 0.499;
    Scalar s_right = 0.501;
    Scalar d2_left = spline.EvaluateSecondDerivative(s_left);
    Scalar d2_right = spline.EvaluateSecondDerivative(s_right);
    EXPECT_NEAR(d2_left, d2_right, 0.01);  // Should be close
}

TEST_F(CubicSplineTest, ExtrapolationBehavior) {
    // Test behavior outside the data range (uses endpoint intervals)
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0};
    std::vector<Scalar> x = {0.0, 1.0, 4.0, 9.0};  // Roughly x = s^2

    CubicSpline spline;
    spline.Build(s, x);

    // Extrapolation below range
    Scalar val_below = spline.Evaluate(-0.5);
    EXPECT_TRUE(std::isfinite(val_below));

    // Extrapolation above range
    Scalar val_above = spline.Evaluate(3.5);
    EXPECT_TRUE(std::isfinite(val_above));
}

TEST_F(CubicSplineTest, InvalidInputThrows) {
    CubicSpline spline;

    // Mismatched sizes should throw
    std::vector<Scalar> s = {0.0, 1.0, 2.0};
    std::vector<Scalar> x = {0.0, 1.0};  // Wrong size
    EXPECT_THROW(spline.Build(s, x), std::invalid_argument);

    // Too few points should throw
    std::vector<Scalar> s_short = {0.0};
    std::vector<Scalar> x_short = {1.0};
    EXPECT_THROW(spline.Build(s_short, x_short), std::invalid_argument);
}

TEST_F(CubicSplineTest, EvaluateBeforeBuildThrows) {
    CubicSpline spline;
    EXPECT_THROW(spline.Evaluate(0.5), std::runtime_error);
    EXPECT_THROW(spline.EvaluateDerivative(0.5), std::runtime_error);
    EXPECT_THROW(spline.EvaluateSecondDerivative(0.5), std::runtime_error);

    Scalar v, d1, d2;
    EXPECT_THROW(spline.EvaluateAll(0.5, v, d1, d2), std::runtime_error);
}

TEST_F(CubicSplineTest, ValidityAndSizeQueries) {
    CubicSpline spline;

    // Before building
    EXPECT_FALSE(spline.IsValid());
    EXPECT_EQ(spline.Size(), 0);

    // After building
    std::vector<Scalar> s = {0.0, 1.0, 2.0, 3.0};
    std::vector<Scalar> x = {1.0, 2.0, 3.0, 4.0};
    spline.Build(s, x);

    EXPECT_TRUE(spline.IsValid());
    EXPECT_EQ(spline.Size(), 4);
    EXPECT_NEAR(spline.MinParameter(), 0.0, kTestTolerance);
    EXPECT_NEAR(spline.MaxParameter(), 3.0, kTestTolerance);
}

}  // namespace
}  // namespace pegasus
