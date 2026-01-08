// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_AD_DUAL_NUMBER_H_
#define VLM_AD_DUAL_NUMBER_H_

#include <cmath>
#include <ostream>
#include <type_traits>

namespace vlm {
namespace ad {

/// Forward-mode automatic differentiation using dual numbers.
///
/// A dual number has the form: a + b*ε where ε² = 0
/// When computing f(a + b*ε), the result is f(a) + f'(a)*b*ε
///
/// This class supports computing derivatives with respect to N variables
/// simultaneously by storing N derivative components.
///
/// @tparam T The underlying scalar type (float, double)
/// @tparam N Number of derivative components (for multi-variable gradients)
///
/// Example usage:
/// @code
///   // Single variable derivative
///   Dual<double, 1> x(2.0, 1.0);  // x = 2, dx/dx = 1
///   auto y = x * x + 3.0 * x;     // y = x² + 3x
///   // y.Value() = 10.0, y.Derivative(0) = 7.0 (dy/dx at x=2)
///
///   // Multi-variable gradient
///   Dual<double, 2> a(1.0, {1.0, 0.0});  // a, ∂/∂a = 1, ∂/∂b = 0
///   Dual<double, 2> b(2.0, {0.0, 1.0});  // b, ∂/∂a = 0, ∂/∂b = 1
///   auto f = a * b + a * a;
///   // f.Value() = 3.0, f.Derivative(0) = 4.0 (∂f/∂a), f.Derivative(1) = 1.0 (∂f/∂b)
/// @endcode
template <typename T = double, int N = 1>
class Dual {
 public:
  static_assert(N >= 1, "Must have at least one derivative component");

  /// Default constructor (zero value and derivatives)
  Dual() : value_(T(0)) { derivatives_.fill(T(0)); }

  /// Construct from value only (derivatives are zero)
  /// This allows implicit conversion from scalars
  Dual(T value) : value_(value) { derivatives_.fill(T(0)); }  // NOLINT

  /// Construct with value and single derivative (for N=1)
  template <int M = N, typename = std::enable_if_t<M == 1>>
  Dual(T value, T derivative) : value_(value) {
    derivatives_[0] = derivative;
  }

  /// Construct with value and derivative array
  Dual(T value, const std::array<T, N>& derivatives)
      : value_(value), derivatives_(derivatives) {}

  /// Construct with value and initializer list for derivatives
  Dual(T value, std::initializer_list<T> derivs) : value_(value) {
    derivatives_.fill(T(0));
    int i = 0;
    for (auto d : derivs) {
      if (i >= N) break;
      derivatives_[i++] = d;
    }
  }

  /// Get the value part
  T Value() const { return value_; }

  /// Get a specific derivative component
  T Derivative(int i = 0) const { return derivatives_[i]; }

  /// Get all derivatives as array
  const std::array<T, N>& Derivatives() const { return derivatives_; }

  /// Set the value
  void SetValue(T v) { value_ = v; }

  /// Set a specific derivative
  void SetDerivative(int i, T d) { derivatives_[i] = d; }

  /// Create a variable for differentiation (sets derivative to 1 for index i)
  static Dual Variable(T value, int derivative_index = 0) {
    Dual result(value);
    if (derivative_index >= 0 && derivative_index < N) {
      result.derivatives_[derivative_index] = T(1);
    }
    return result;
  }

  /// Create a constant (all derivatives zero)
  static Dual Constant(T value) { return Dual(value); }

  // Arithmetic operators

  Dual operator+() const { return *this; }

  Dual operator-() const {
    std::array<T, N> neg_derivs;
    for (int i = 0; i < N; ++i) {
      neg_derivs[i] = -derivatives_[i];
    }
    return Dual(-value_, neg_derivs);
  }

  Dual& operator+=(const Dual& other) {
    value_ += other.value_;
    for (int i = 0; i < N; ++i) {
      derivatives_[i] += other.derivatives_[i];
    }
    return *this;
  }

  Dual& operator-=(const Dual& other) {
    value_ -= other.value_;
    for (int i = 0; i < N; ++i) {
      derivatives_[i] -= other.derivatives_[i];
    }
    return *this;
  }

  Dual& operator*=(const Dual& other) {
    // d(a*b) = a*db + b*da
    for (int i = 0; i < N; ++i) {
      derivatives_[i] =
          value_ * other.derivatives_[i] + derivatives_[i] * other.value_;
    }
    value_ *= other.value_;
    return *this;
  }

  Dual& operator/=(const Dual& other) {
    // d(a/b) = (da*b - a*db) / b²
    T inv_b = T(1) / other.value_;
    T inv_b2 = inv_b * inv_b;
    for (int i = 0; i < N; ++i) {
      derivatives_[i] =
          (derivatives_[i] * other.value_ - value_ * other.derivatives_[i]) *
          inv_b2;
    }
    value_ *= inv_b;
    return *this;
  }

  // Scalar arithmetic (right-hand side)
  Dual& operator+=(T s) {
    value_ += s;
    return *this;
  }

  Dual& operator-=(T s) {
    value_ -= s;
    return *this;
  }

  Dual& operator*=(T s) {
    value_ *= s;
    for (int i = 0; i < N; ++i) {
      derivatives_[i] *= s;
    }
    return *this;
  }

  Dual& operator/=(T s) {
    T inv = T(1) / s;
    value_ *= inv;
    for (int i = 0; i < N; ++i) {
      derivatives_[i] *= inv;
    }
    return *this;
  }

  // Comparison operators (compare values only)
  bool operator==(const Dual& other) const { return value_ == other.value_; }
  bool operator!=(const Dual& other) const { return value_ != other.value_; }
  bool operator<(const Dual& other) const { return value_ < other.value_; }
  bool operator<=(const Dual& other) const { return value_ <= other.value_; }
  bool operator>(const Dual& other) const { return value_ > other.value_; }
  bool operator>=(const Dual& other) const { return value_ >= other.value_; }

  // Comparison with scalars
  bool operator==(T s) const { return value_ == s; }
  bool operator!=(T s) const { return value_ != s; }
  bool operator<(T s) const { return value_ < s; }
  bool operator<=(T s) const { return value_ <= s; }
  bool operator>(T s) const { return value_ > s; }
  bool operator>=(T s) const { return value_ >= s; }

  // Conversion to underlying type (returns value only)
  explicit operator T() const { return value_; }

 private:
  T value_;
  std::array<T, N> derivatives_;
};

// Binary operators

template <typename T, int N>
Dual<T, N> operator+(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> result = a;
  result += b;
  return result;
}

template <typename T, int N>
Dual<T, N> operator-(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> result = a;
  result -= b;
  return result;
}

template <typename T, int N>
Dual<T, N> operator*(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> result = a;
  result *= b;
  return result;
}

template <typename T, int N>
Dual<T, N> operator/(const Dual<T, N>& a, const Dual<T, N>& b) {
  Dual<T, N> result = a;
  result /= b;
  return result;
}

// Scalar-dual operations (left-hand scalar)

template <typename T, int N>
Dual<T, N> operator+(T s, const Dual<T, N>& d) {
  return d + s;
}

template <typename T, int N>
Dual<T, N> operator-(T s, const Dual<T, N>& d) {
  return Dual<T, N>(s) - d;
}

template <typename T, int N>
Dual<T, N> operator*(T s, const Dual<T, N>& d) {
  Dual<T, N> result = d;
  result *= s;
  return result;
}

template <typename T, int N>
Dual<T, N> operator/(T s, const Dual<T, N>& d) {
  return Dual<T, N>(s) / d;
}

// Scalar-dual operations (right-hand scalar)

template <typename T, int N>
Dual<T, N> operator+(const Dual<T, N>& d, T s) {
  Dual<T, N> result = d;
  result += s;
  return result;
}

template <typename T, int N>
Dual<T, N> operator-(const Dual<T, N>& d, T s) {
  Dual<T, N> result = d;
  result -= s;
  return result;
}

template <typename T, int N>
Dual<T, N> operator*(const Dual<T, N>& d, T s) {
  Dual<T, N> result = d;
  result *= s;
  return result;
}

template <typename T, int N>
Dual<T, N> operator/(const Dual<T, N>& d, T s) {
  Dual<T, N> result = d;
  result /= s;
  return result;
}

// Mathematical functions

template <typename T, int N>
Dual<T, N> abs(const Dual<T, N>& x) {
  if (x.Value() >= T(0)) {
    return x;
  } else {
    return -x;
  }
}

template <typename T, int N>
Dual<T, N> sqrt(const Dual<T, N>& x) {
  T val = std::sqrt(x.Value());
  T deriv_coeff = T(0.5) / val;  // d/dx sqrt(x) = 1/(2*sqrt(x))
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> sin(const Dual<T, N>& x) {
  T val = std::sin(x.Value());
  T deriv_coeff = std::cos(x.Value());  // d/dx sin(x) = cos(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> cos(const Dual<T, N>& x) {
  T val = std::cos(x.Value());
  T deriv_coeff = -std::sin(x.Value());  // d/dx cos(x) = -sin(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> tan(const Dual<T, N>& x) {
  T val = std::tan(x.Value());
  T c = std::cos(x.Value());
  T deriv_coeff = T(1) / (c * c);  // d/dx tan(x) = sec²(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> asin(const Dual<T, N>& x) {
  T val = std::asin(x.Value());
  // d/dx asin(x) = 1/sqrt(1-x²)
  T deriv_coeff = T(1) / std::sqrt(T(1) - x.Value() * x.Value());
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> acos(const Dual<T, N>& x) {
  T val = std::acos(x.Value());
  // d/dx acos(x) = -1/sqrt(1-x²)
  T deriv_coeff = T(-1) / std::sqrt(T(1) - x.Value() * x.Value());
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> atan(const Dual<T, N>& x) {
  T val = std::atan(x.Value());
  // d/dx atan(x) = 1/(1+x²)
  T deriv_coeff = T(1) / (T(1) + x.Value() * x.Value());
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> atan2(const Dual<T, N>& y, const Dual<T, N>& x) {
  T val = std::atan2(y.Value(), x.Value());
  // d/dy atan2(y,x) = x/(x²+y²), d/dx atan2(y,x) = -y/(x²+y²)
  T r2 = x.Value() * x.Value() + y.Value() * y.Value();
  T dx_coeff = -y.Value() / r2;
  T dy_coeff = x.Value() / r2;
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = dy_coeff * y.Derivative(i) + dx_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> exp(const Dual<T, N>& x) {
  T val = std::exp(x.Value());
  // d/dx exp(x) = exp(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = val * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> log(const Dual<T, N>& x) {
  T val = std::log(x.Value());
  // d/dx log(x) = 1/x
  T deriv_coeff = T(1) / x.Value();
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> log10(const Dual<T, N>& x) {
  T val = std::log10(x.Value());
  // d/dx log10(x) = 1/(x*ln(10))
  T deriv_coeff = T(1) / (x.Value() * std::log(T(10)));
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> pow(const Dual<T, N>& base, const Dual<T, N>& exponent) {
  T val = std::pow(base.Value(), exponent.Value());
  // d/dx a^b = a^b * (b*da/a + log(a)*db)
  T log_base = std::log(base.Value());
  T b_over_a = exponent.Value() / base.Value();
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = val * (b_over_a * base.Derivative(i) +
                       log_base * exponent.Derivative(i));
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> pow(const Dual<T, N>& base, T exponent) {
  T val = std::pow(base.Value(), exponent);
  // d/dx a^n = n * a^(n-1) * da
  T deriv_coeff = exponent * std::pow(base.Value(), exponent - T(1));
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * base.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> pow(T base, const Dual<T, N>& exponent) {
  T val = std::pow(base, exponent.Value());
  // d/dx c^b = c^b * log(c) * db
  T deriv_coeff = val * std::log(base);
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * exponent.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> sinh(const Dual<T, N>& x) {
  T val = std::sinh(x.Value());
  T deriv_coeff = std::cosh(x.Value());  // d/dx sinh(x) = cosh(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> cosh(const Dual<T, N>& x) {
  T val = std::cosh(x.Value());
  T deriv_coeff = std::sinh(x.Value());  // d/dx cosh(x) = sinh(x)
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> tanh(const Dual<T, N>& x) {
  T val = std::tanh(x.Value());
  // d/dx tanh(x) = 1 - tanh²(x) = sech²(x)
  T deriv_coeff = T(1) - val * val;
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = deriv_coeff * x.Derivative(i);
  }
  return Dual<T, N>(val, derivs);
}

template <typename T, int N>
Dual<T, N> hypot(const Dual<T, N>& x, const Dual<T, N>& y) {
  T val = std::hypot(x.Value(), y.Value());
  // d/dx sqrt(x²+y²) = x/sqrt(x²+y²), similarly for y
  T inv_val = T(1) / val;
  std::array<T, N> derivs;
  for (int i = 0; i < N; ++i) {
    derivs[i] = inv_val * (x.Value() * x.Derivative(i) +
                           y.Value() * y.Derivative(i));
  }
  return Dual<T, N>(val, derivs);
}

// Min/max (discontinuous derivatives - uses subgradient)

template <typename T, int N>
Dual<T, N> fmin(const Dual<T, N>& a, const Dual<T, N>& b) {
  return (a.Value() <= b.Value()) ? a : b;
}

template <typename T, int N>
Dual<T, N> fmax(const Dual<T, N>& a, const Dual<T, N>& b) {
  return (a.Value() >= b.Value()) ? a : b;
}

// Floor/ceil (derivative is zero almost everywhere)

template <typename T, int N>
Dual<T, N> floor(const Dual<T, N>& x) {
  return Dual<T, N>(std::floor(x.Value()));
}

template <typename T, int N>
Dual<T, N> ceil(const Dual<T, N>& x) {
  return Dual<T, N>(std::ceil(x.Value()));
}

// Stream output

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, const Dual<T, N>& d) {
  os << d.Value();
  if constexpr (N == 1) {
    os << " + " << d.Derivative(0) << "ε";
  } else {
    os << " + [";
    for (int i = 0; i < N; ++i) {
      if (i > 0) os << ", ";
      os << d.Derivative(i);
    }
    os << "]ε";
  }
  return os;
}

// Convenience type aliases

/// Single derivative component (for computing df/dx)
using Dual1d = Dual<double, 1>;

/// Two derivative components (for gradient of function of 2 variables)
using Dual2d = Dual<double, 2>;

/// Three derivative components (for gradient of function of 3 variables)
using Dual3d = Dual<double, 3>;

/// Four derivative components
using Dual4d = Dual<double, 4>;

/// Six derivative components (common for 6 DOF systems)
using Dual6d = Dual<double, 6>;

}  // namespace ad
}  // namespace vlm

// Eigen NumTraits specialization for use in Eigen matrices
// This allows Dual numbers to be used in Eigen::Matrix<Dual<T,N>, ...>
#include <Eigen/Core>

namespace Eigen {

template <typename T, int N>
struct NumTraits<vlm::ad::Dual<T, N>> : NumTraits<T> {
  using Real = vlm::ad::Dual<T, N>;
  using NonInteger = vlm::ad::Dual<T, N>;
  using Nested = vlm::ad::Dual<T, N>;
  using Literal = T;

  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = N + 1,
    AddCost = N + 1,
    MulCost = 2 * N + 1
  };

  static inline Real epsilon() {
    return Real(NumTraits<T>::epsilon());
  }
  static inline Real dummy_precision() {
    return Real(NumTraits<T>::dummy_precision());
  }
  static inline int digits10() { return NumTraits<T>::digits10(); }
};

}  // namespace Eigen

#endif  // VLM_AD_DUAL_NUMBER_H_
