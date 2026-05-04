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

#ifndef PEGASUS_UTILS_CUBIC_SPLINE_H_
#define PEGASUS_UTILS_CUBIC_SPLINE_H_

#include "pegasus/core/types.h"

#include <vector>

namespace pegasus {

/// @brief Boundary condition type for cubic spline endpoints.
enum class SplineBoundaryCondition {
    kNatural,         ///< Zero second derivative (natural spline)
    kZeroThirdDeriv,  ///< Zero third derivative
    kSpecifiedSlope   ///< Specified first derivative value
};

/// @brief Cubic spline interpolator for 1D data.
///
/// Implements cubic spline interpolation with configurable boundary conditions.
/// The spline is parameterized by an independent variable s, and interpolates
/// dependent variable x(s).
///
/// Uses the algorithm from AVL's spline.f, which solves a tridiagonal system
/// to compute spline slopes at each data point.
///
/// @par Example usage:
/// @code
/// std::vector<double> s = {0.0, 1.0, 2.0, 3.0};
/// std::vector<double> x = {0.0, 1.0, 0.0, 1.0};
/// CubicSpline spline;
/// spline.Build(s, x);
/// double value = spline.Evaluate(1.5);
/// double slope = spline.EvaluateDerivative(1.5);
/// @endcode
///
/// Reference: Drela, M. & Youngren, H., AVL spline.f
class CubicSpline {
public:
    CubicSpline() = default;

    /// @brief Builds the spline from data with natural boundary conditions.
    ///
    /// Natural boundary conditions impose zero second derivative at both ends.
    ///
    /// @param s Independent variable values (must be monotonically increasing)
    /// @param x Dependent variable values
    /// @throws std::invalid_argument if s and x have different sizes or size <
    /// 2
    void Build(const std::vector<Scalar>& s, const std::vector<Scalar>& x);

    /// @brief Builds the spline with specified boundary conditions.
    ///
    /// @param s Independent variable values (must be monotonically increasing)
    /// @param x Dependent variable values
    /// @param left_condition Boundary condition at left endpoint
    /// @param right_condition Boundary condition at right endpoint
    /// @param left_slope First derivative at left endpoint (used if
    ///        left_condition is kSpecifiedSlope)
    /// @param right_slope First derivative at right endpoint (used if
    ///        right_condition is kSpecifiedSlope)
    void Build(const std::vector<Scalar>& s, const std::vector<Scalar>& x,
               SplineBoundaryCondition left_condition,
               SplineBoundaryCondition right_condition, Scalar left_slope = 0.0,
               Scalar right_slope = 0.0);

    /// @brief Evaluates the spline at a given parameter value.
    ///
    /// @param s_eval Parameter value at which to evaluate
    /// @return Interpolated value x(s_eval)
    [[nodiscard]] Scalar Evaluate(Scalar s_eval) const;

    /// @brief Evaluates the first derivative dx/ds at a given parameter value.
    ///
    /// @param s_eval Parameter value at which to evaluate
    /// @return First derivative dx/ds at s_eval
    [[nodiscard]] Scalar EvaluateDerivative(Scalar s_eval) const;

    /// @brief Evaluates the second derivative d^2x/ds^2 at a given parameter.
    ///
    /// @param s_eval Parameter value at which to evaluate
    /// @return Second derivative d^2x/ds^2 at s_eval
    [[nodiscard]] Scalar EvaluateSecondDerivative(Scalar s_eval) const;

    /// @brief Evaluates value and all derivatives at once (more efficient).
    ///
    /// @param s_eval Parameter value at which to evaluate
    /// @param[out] value Interpolated value x(s_eval)
    /// @param[out] first_deriv First derivative dx/ds
    /// @param[out] second_deriv Second derivative d^2x/ds^2
    void EvaluateAll(Scalar s_eval, Scalar& value, Scalar& first_deriv,
                     Scalar& second_deriv) const;

    /// @brief Returns the number of data points.
    [[nodiscard]] Index Size() const { return static_cast<Index>(s_.size()); }

    /// @brief Checks if the spline has been built with data.
    [[nodiscard]] bool IsValid() const { return !s_.empty(); }

    /// @brief Returns the minimum parameter value.
    [[nodiscard]] Scalar MinParameter() const {
        return s_.empty() ? 0.0 : s_.front();
    }

    /// @brief Returns the maximum parameter value.
    [[nodiscard]] Scalar MaxParameter() const {
        return s_.empty() ? 0.0 : s_.back();
    }

private:
    /// @brief Finds the interval index for a given parameter value.
    ///
    /// Uses binary search to find i such that s_[i-1] <= s_eval < s_[i].
    ///
    /// @param s_eval Parameter value to locate
    /// @return Index i of the upper bound of the containing interval
    [[nodiscard]] Index FindInterval(Scalar s_eval) const;

    /// @brief Solves the tridiagonal system for spline slopes.
    ///
    /// Solves the system:
    ///   | a[0] c[0]                   | | xs[0]   |   | d[0]   |
    ///   | b[1] a[1] c[1]              | | xs[1]   |   | d[1]   |
    ///   |      b[2] a[2] c[2]         | | xs[2]   | = | d[2]   |
    ///   |           ...               | | ...     |   | ...    |
    ///   |                b[n-1] a[n-1]| | xs[n-1] |   | d[n-1] |
    ///
    /// @param a Diagonal coefficients (destroyed on output)
    /// @param b Sub-diagonal coefficients
    /// @param c Super-diagonal coefficients (destroyed on output)
    /// @param d Right-hand side (replaced by solution on output)
    static void SolveTridiagonal(std::vector<Scalar>& a,
                                 const std::vector<Scalar>& b,
                                 std::vector<Scalar>& c,
                                 std::vector<Scalar>& d);

    std::vector<Scalar> s_;   ///< Independent variable values
    std::vector<Scalar> x_;   ///< Dependent variable values
    std::vector<Scalar> xs_;  ///< Spline slopes dx/ds at each point
};

}  // namespace pegasus

#endif  // PEGASUS_UTILS_CUBIC_SPLINE_H_
