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

#include <algorithm>
#include <stdexcept>

namespace pegasus {

void CubicSpline::Build(const std::vector<Scalar>& s,
                        const std::vector<Scalar>& x) {
    Build(s, x, SplineBoundaryCondition::kNatural,
          SplineBoundaryCondition::kNatural);
}

void CubicSpline::Build(const std::vector<Scalar>& s,
                        const std::vector<Scalar>& x,
                        SplineBoundaryCondition left_condition,
                        SplineBoundaryCondition right_condition,
                        Scalar left_slope, Scalar right_slope) {
    if (s.size() != x.size()) {
        throw std::invalid_argument(
            "CubicSpline::Build: s and x must have the same size");
    }
    if (s.size() < 2) {
        throw std::invalid_argument(
            "CubicSpline::Build: need at least 2 data points");
    }

    const auto n = static_cast<Index>(s.size());
    s_ = s;
    x_ = x;
    xs_.resize(n);

    // Working arrays for tridiagonal system
    std::vector<Scalar> a(n);  // Diagonal
    std::vector<Scalar> b(n);  // Sub-diagonal
    std::vector<Scalar> c(n);  // Super-diagonal

    // Set up interior equations
    // For i = 1 to n-2 (0-indexed), the continuity conditions give:
    // ds_minus * xs[i-1] + 2*(ds_minus + ds_plus) * xs[i] + ds_plus * xs[i+1]
    //   = 3 * ((x[i+1]-x[i])*ds_minus/ds_plus + (x[i]-x[i-1])*ds_plus/ds_minus)
    for (Index i = 1; i < n - 1; ++i) {
        Scalar ds_minus = s_[i] - s_[i - 1];
        Scalar ds_plus = s_[i + 1] - s_[i];
        b[i] = ds_plus;
        a[i] = 2.0 * (ds_minus + ds_plus);
        c[i] = ds_minus;
        xs_[i] = 3.0 * ((x_[i + 1] - x_[i]) * ds_minus / ds_plus +
                        (x_[i] - x_[i - 1]) * ds_plus / ds_minus);
    }

    // Set left boundary condition
    Scalar ds_left = s_[1] - s_[0];
    switch (left_condition) {
        case SplineBoundaryCondition::kNatural:
            // Zero second derivative: xs[0] + 0.5*xs[1] = 1.5*(x[1]-x[0])/ds
            a[0] = 2.0;
            c[0] = 1.0;
            xs_[0] = 3.0 * (x_[1] - x_[0]) / ds_left;
            break;
        case SplineBoundaryCondition::kZeroThirdDeriv:
            // Zero third derivative
            a[0] = 1.0;
            c[0] = 1.0;
            xs_[0] = 2.0 * (x_[1] - x_[0]) / ds_left;
            break;
        case SplineBoundaryCondition::kSpecifiedSlope:
            // Specified first derivative
            a[0] = 1.0;
            c[0] = 0.0;
            xs_[0] = left_slope;
            break;
    }

    // Set right boundary condition
    Scalar ds_right = s_[n - 1] - s_[n - 2];
    switch (right_condition) {
        case SplineBoundaryCondition::kNatural:
            // Zero second derivative
            b[n - 1] = 1.0;
            a[n - 1] = 2.0;
            xs_[n - 1] = 3.0 * (x_[n - 1] - x_[n - 2]) / ds_right;
            break;
        case SplineBoundaryCondition::kZeroThirdDeriv:
            // Zero third derivative
            b[n - 1] = 1.0;
            a[n - 1] = 1.0;
            xs_[n - 1] = 2.0 * (x_[n - 1] - x_[n - 2]) / ds_right;
            break;
        case SplineBoundaryCondition::kSpecifiedSlope:
            // Specified first derivative
            a[n - 1] = 1.0;
            b[n - 1] = 0.0;
            xs_[n - 1] = right_slope;
            break;
    }

    // Special case: only 2 points with zero third derivative at both ends
    // is not well-posed, fall back to natural spline at right end
    if (n == 2 && left_condition == SplineBoundaryCondition::kZeroThirdDeriv &&
        right_condition == SplineBoundaryCondition::kZeroThirdDeriv) {
        b[n - 1] = 1.0;
        a[n - 1] = 2.0;
        xs_[n - 1] = 3.0 * (x_[n - 1] - x_[n - 2]) / ds_right;
    }

    // Solve the tridiagonal system for slopes xs_
    SolveTridiagonal(a, b, c, xs_);
}

Index CubicSpline::FindInterval(Scalar s_eval) const {
    // Binary search to find interval
    // Returns i such that s_[i-1] <= s_eval < s_[i]
    // For values outside the range, clamp to valid interval

    const auto n = static_cast<Index>(s_.size());

    if (s_eval <= s_[0]) {
        return 1;
    }
    if (s_eval >= s_[n - 1]) {
        return n - 1;
    }

    Index low = 0;
    Index high = n - 1;

    while (high - low > 1) {
        Index mid = (low + high) / 2;
        if (s_eval < s_[mid]) {
            high = mid;
        } else {
            low = mid;
        }
    }

    return high;
}

Scalar CubicSpline::Evaluate(Scalar s_eval) const {
    if (!IsValid()) {
        throw std::runtime_error("CubicSpline::Evaluate: spline not built");
    }

    Index i = FindInterval(s_eval);
    Scalar ds = s_[i] - s_[i - 1];
    Scalar t = (s_eval - s_[i - 1]) / ds;

    // Cubic Hermite interpolation
    // x(t) = (1-t)*x[i-1] + t*x[i] + t*(1-t)*[(1-t)*cx1 - t*cx2]
    // where:
    //   cx1 = ds*xs[i-1] - (x[i] - x[i-1])
    //   cx2 = ds*xs[i]   - (x[i] - x[i-1])
    Scalar dx = x_[i] - x_[i - 1];
    Scalar cx1 = ds * xs_[i - 1] - dx;
    Scalar cx2 = ds * xs_[i] - dx;

    return t * x_[i] + (1.0 - t) * x_[i - 1] +
           t * (1.0 - t) * ((1.0 - t) * cx1 - t * cx2);
}

Scalar CubicSpline::EvaluateDerivative(Scalar s_eval) const {
    if (!IsValid()) {
        throw std::runtime_error(
            "CubicSpline::EvaluateDerivative: spline not built");
    }

    Index i = FindInterval(s_eval);
    Scalar ds = s_[i] - s_[i - 1];
    Scalar t = (s_eval - s_[i - 1]) / ds;

    Scalar dx = x_[i] - x_[i - 1];
    Scalar cx1 = ds * xs_[i - 1] - dx;
    Scalar cx2 = ds * xs_[i] - dx;

    // Derivative of cubic Hermite formula with respect to t, then divide by ds
    Scalar dxdt =
        dx + (1.0 - 4.0 * t + 3.0 * t * t) * cx1 + t * (3.0 * t - 2.0) * cx2;

    return dxdt / ds;
}

Scalar CubicSpline::EvaluateSecondDerivative(Scalar s_eval) const {
    if (!IsValid()) {
        throw std::runtime_error(
            "CubicSpline::EvaluateSecondDerivative: spline not built");
    }

    Index i = FindInterval(s_eval);
    Scalar ds = s_[i] - s_[i - 1];
    Scalar t = (s_eval - s_[i - 1]) / ds;

    Scalar dx = x_[i] - x_[i - 1];
    Scalar cx1 = ds * xs_[i - 1] - dx;
    Scalar cx2 = ds * xs_[i] - dx;

    // Second derivative with respect to t, then divide by ds^2
    Scalar d2xdt2 = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

    return d2xdt2 / (ds * ds);
}

void CubicSpline::EvaluateAll(Scalar s_eval, Scalar& value, Scalar& first_deriv,
                              Scalar& second_deriv) const {
    if (!IsValid()) {
        throw std::runtime_error("CubicSpline::EvaluateAll: spline not built");
    }

    Index i = FindInterval(s_eval);
    Scalar ds = s_[i] - s_[i - 1];
    Scalar t = (s_eval - s_[i - 1]) / ds;

    // Alternative form using polynomial coefficients (more efficient)
    // x(t) = f0 + t*f1 + t^2*f2 + t^3*f3
    Scalar f0 = x_[i - 1];
    Scalar f1 = ds * xs_[i - 1];
    Scalar f2 = -ds * (2.0 * xs_[i - 1] + xs_[i]) + 3.0 * (x_[i] - x_[i - 1]);
    Scalar f3 = ds * (xs_[i - 1] + xs_[i]) - 2.0 * (x_[i] - x_[i - 1]);

    value = f0 + t * (f1 + t * (f2 + t * f3));
    first_deriv = (f1 + t * (2.0 * f2 + t * 3.0 * f3)) / ds;
    second_deriv = (2.0 * f2 + t * 6.0 * f3) / (ds * ds);
}

void CubicSpline::SolveTridiagonal(std::vector<Scalar>& a,
                                   const std::vector<Scalar>& b,
                                   std::vector<Scalar>& c,
                                   std::vector<Scalar>& d) {
    const auto n = static_cast<Index>(a.size());

    // Forward elimination
    for (Index k = 1; k < n; ++k) {
        c[k - 1] = c[k - 1] / a[k - 1];
        d[k - 1] = d[k - 1] / a[k - 1];
        a[k] = a[k] - b[k] * c[k - 1];
        d[k] = d[k] - b[k] * d[k - 1];
    }

    // Last element
    d[n - 1] = d[n - 1] / a[n - 1];

    // Back substitution
    for (Index k = n - 2; k >= 0; --k) {
        d[k] = d[k] - c[k] * d[k + 1];
    }
}

}  // namespace pegasus
