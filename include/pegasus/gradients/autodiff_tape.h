/*
    Project Pegasus
    Version 0.1.0

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

#ifndef PEGASUS_GRADIENTS_AUTODIFF_TAPE_H_
#define PEGASUS_GRADIENTS_AUTODIFF_TAPE_H_

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace pegasus {
namespace gradients {

/// Forward declaration
class TapeVariable;
class Tape;

/// Global tape for recording operations (thread-local for safety)
inline thread_local Tape* active_tape = nullptr;

/// Operation codes for the tape
enum class OpCode {
    kConstant,  // No operation, just a constant
    kAdd,       // a + b
    kSub,       // a - b
    kMul,       // a * b
    kDiv,       // a / b
    kNeg,       // -a
    kSqrt,      // sqrt(a)
    kSin,       // sin(a)
    kCos,       // cos(a)
    kTan,       // tan(a)
    kExp,       // exp(a)
    kLog,       // log(a)
    kPow,       // pow(a, b)
    kAbs,       // abs(a)
    kAsin,      // asin(a)
    kAcos,      // acos(a)
    kAtan,      // atan(a)
    kAtan2,     // atan2(a, b)
    kSinh,      // sinh(a)
    kCosh,      // cosh(a)
    kTanh,      // tanh(a)
};

/// Entry in the computation tape
struct TapeEntry {
    OpCode op;
    double value;     // Result of the operation
    int operand1;     // Index of first operand (-1 if none)
    int operand2;     // Index of second operand (-1 if none)
    double partial1;  // Partial derivative w.r.t. operand1
    double partial2;  // Partial derivative w.r.t. operand2
};

/// Computation tape for reverse-mode AD
///
/// Records all operations during forward evaluation, then computes
/// gradients by traversing the tape backwards.
///
/// Example usage:
/// @code
///   Tape tape;
///   tape.Activate();
///
///   TapeVariable x = tape.CreateVariable(2.0);
///   TapeVariable y = tape.CreateVariable(3.0);
///   TapeVariable z = x * y + sin(x);
///
///   auto gradients = tape.ComputeGradient(z);
///   // gradients[x.Index()] = dz/dx = y + cos(x)
///   // gradients[y.Index()] = dz/dy = x
///
///   tape.Deactivate();
/// @endcode
class Tape {
public:
    Tape() = default;

    /// Activate this tape for recording
    void Activate() { active_tape = this; }

    /// Deactivate this tape
    void Deactivate() {
        if (active_tape == this) {
            active_tape = nullptr;
        }
    }

    /// Clear all recorded operations
    void Clear() {
        entries_.clear();
        num_variables_ = 0;
    }

    /// Create an input variable
    TapeVariable CreateVariable(double value);

    /// Create a constant (not differentiated)
    TapeVariable CreateConstant(double value);

    /// Record an operation on the tape
    int RecordOperation(OpCode op, double value, int op1, int op2, double p1,
                        double p2) {
        int idx = static_cast<int>(entries_.size());
        entries_.push_back({op, value, op1, op2, p1, p2});
        return idx;
    }

    /// Get entry by index
    const TapeEntry& GetEntry(int index) const { return entries_[index]; }

    /// Get number of entries
    int NumEntries() const { return static_cast<int>(entries_.size()); }

    /// Get number of input variables
    int NumVariables() const { return num_variables_; }

    /// Compute gradient of output w.r.t. all variables
    std::vector<double> ComputeGradient(const TapeVariable& output) const;

    /// Compute gradient of output w.r.t. specified variable
    double ComputeGradient(const TapeVariable& output,
                           const TapeVariable& variable) const;

private:
    std::vector<TapeEntry> entries_;
    int num_variables_ = 0;

    friend class TapeVariable;
};

/// Variable that records operations to a tape
///
/// This class represents a value in a differentiable computation.
/// Operations on TapeVariable instances are recorded to the active tape.
class TapeVariable {
public:
    /// Default constructor (creates uninitialized variable)
    TapeVariable() : tape_(nullptr), index_(-1), value_(0.0) {}

    /// Get the value
    double Value() const { return value_; }

    /// Get the tape index
    int Index() const { return index_; }

    /// Check if valid
    bool IsValid() const { return tape_ != nullptr && index_ >= 0; }

    /// Get the tape this variable belongs to
    Tape* GetTape() const { return tape_; }

    // Arithmetic operators

    TapeVariable operator+() const { return *this; }

    TapeVariable operator-() const {
        if (!tape_)
            return TapeVariable(nullptr, -1, -value_);
        int idx = tape_->RecordOperation(OpCode::kNeg, -value_, index_, -1,
                                         -1.0, 0.0);
        return TapeVariable(tape_, idx, -value_);
    }

    TapeVariable& operator+=(const TapeVariable& other);
    TapeVariable& operator-=(const TapeVariable& other);
    TapeVariable& operator*=(const TapeVariable& other);
    TapeVariable& operator/=(const TapeVariable& other);

    TapeVariable& operator+=(double s);
    TapeVariable& operator-=(double s);
    TapeVariable& operator*=(double s);
    TapeVariable& operator/=(double s);

    // Comparison operators (compare values only)
    bool operator==(const TapeVariable& other) const {
        return value_ == other.value_;
    }
    bool operator!=(const TapeVariable& other) const {
        return value_ != other.value_;
    }
    bool operator<(const TapeVariable& other) const {
        return value_ < other.value_;
    }
    bool operator<=(const TapeVariable& other) const {
        return value_ <= other.value_;
    }
    bool operator>(const TapeVariable& other) const {
        return value_ > other.value_;
    }
    bool operator>=(const TapeVariable& other) const {
        return value_ >= other.value_;
    }

    bool operator==(double s) const { return value_ == s; }
    bool operator!=(double s) const { return value_ != s; }
    bool operator<(double s) const { return value_ < s; }
    bool operator<=(double s) const { return value_ <= s; }
    bool operator>(double s) const { return value_ > s; }
    bool operator>=(double s) const { return value_ >= s; }

    // Conversion to double
    explicit operator double() const { return value_; }

private:
    TapeVariable(Tape* tape, int index, double value)
        : tape_(tape), index_(index), value_(value) {}

    Tape* tape_;
    int index_;
    double value_;

    friend class Tape;

    // Binary operator friends
    friend TapeVariable operator+(const TapeVariable& a, const TapeVariable& b);
    friend TapeVariable operator-(const TapeVariable& a, const TapeVariable& b);
    friend TapeVariable operator*(const TapeVariable& a, const TapeVariable& b);
    friend TapeVariable operator/(const TapeVariable& a, const TapeVariable& b);

    // Scalar-variable operator friends
    friend TapeVariable operator+(double s, const TapeVariable& v);
    friend TapeVariable operator-(double s, const TapeVariable& v);
    friend TapeVariable operator*(double s, const TapeVariable& v);
    friend TapeVariable operator/(double s, const TapeVariable& v);
    friend TapeVariable operator+(const TapeVariable& v, double s);
    friend TapeVariable operator-(const TapeVariable& v, double s);
    friend TapeVariable operator*(const TapeVariable& v, double s);
    friend TapeVariable operator/(const TapeVariable& v, double s);

    // Math function friends
    friend TapeVariable sqrt(const TapeVariable& x);
    friend TapeVariable sin(const TapeVariable& x);
    friend TapeVariable cos(const TapeVariable& x);
    friend TapeVariable tan(const TapeVariable& x);
    friend TapeVariable exp(const TapeVariable& x);
    friend TapeVariable log(const TapeVariable& x);
    friend TapeVariable pow(const TapeVariable& base, const TapeVariable& exp);
    friend TapeVariable pow(const TapeVariable& base, double exp);
    friend TapeVariable pow(double base, const TapeVariable& exp);
    friend TapeVariable abs(const TapeVariable& x);
    friend TapeVariable asin(const TapeVariable& x);
    friend TapeVariable acos(const TapeVariable& x);
    friend TapeVariable atan(const TapeVariable& x);
    friend TapeVariable atan2(const TapeVariable& y, const TapeVariable& x);
    friend TapeVariable sinh(const TapeVariable& x);
    friend TapeVariable cosh(const TapeVariable& x);
    friend TapeVariable tanh(const TapeVariable& x);
};

// Tape implementation

inline TapeVariable Tape::CreateVariable(double value) {
    int idx = static_cast<int>(entries_.size());
    entries_.push_back({OpCode::kConstant, value, -1, -1, 0.0, 0.0});
    ++num_variables_;
    return TapeVariable(this, idx, value);
}

inline TapeVariable Tape::CreateConstant(double value) {
    int idx = static_cast<int>(entries_.size());
    entries_.push_back({OpCode::kConstant, value, -1, -1, 0.0, 0.0});
    return TapeVariable(this, idx, value);
}

inline std::vector<double> Tape::ComputeGradient(
    const TapeVariable& output) const {
    int n = static_cast<int>(entries_.size());
    std::vector<double> adjoints(n, 0.0);

    // Seed the output adjoint
    if (output.index_ >= 0 && output.index_ < n) {
        adjoints[output.index_] = 1.0;
    }

    // Traverse tape backwards
    for (int i = n - 1; i >= 0; --i) {
        const TapeEntry& entry = entries_[i];
        double adj = adjoints[i];

        if (adj == 0.0)
            continue;

        // Propagate adjoint to operands
        if (entry.operand1 >= 0) {
            adjoints[entry.operand1] += adj * entry.partial1;
        }
        if (entry.operand2 >= 0) {
            adjoints[entry.operand2] += adj * entry.partial2;
        }
    }

    // Extract gradients for input variables (first num_variables_ entries)
    return std::vector<double>(adjoints.begin(),
                               adjoints.begin() + num_variables_);
}

inline double Tape::ComputeGradient(const TapeVariable& output,
                                    const TapeVariable& variable) const {
    auto gradients = ComputeGradient(output);
    if (variable.index_ >= 0 &&
        variable.index_ < static_cast<int>(gradients.size())) {
        return gradients[variable.index_];
    }
    return 0.0;
}

// TapeVariable compound assignment operators

inline TapeVariable& TapeVariable::operator+=(const TapeVariable& other) {
    *this = *this + other;
    return *this;
}

inline TapeVariable& TapeVariable::operator-=(const TapeVariable& other) {
    *this = *this - other;
    return *this;
}

inline TapeVariable& TapeVariable::operator*=(const TapeVariable& other) {
    *this = *this * other;
    return *this;
}

inline TapeVariable& TapeVariable::operator/=(const TapeVariable& other) {
    *this = *this / other;
    return *this;
}

inline TapeVariable& TapeVariable::operator+=(double s) {
    *this = *this + s;
    return *this;
}

inline TapeVariable& TapeVariable::operator-=(double s) {
    *this = *this - s;
    return *this;
}

inline TapeVariable& TapeVariable::operator*=(double s) {
    *this = *this * s;
    return *this;
}

inline TapeVariable& TapeVariable::operator/=(double s) {
    *this = *this / s;
    return *this;
}

// Binary operators

inline TapeVariable operator+(const TapeVariable& a, const TapeVariable& b) {
    Tape* tape = a.tape_ ? a.tape_ : b.tape_;
    if (!tape)
        return TapeVariable(nullptr, -1, a.value_ + b.value_);

    double val = a.value_ + b.value_;
    int idx =
        tape->RecordOperation(OpCode::kAdd, val, a.index_, b.index_, 1.0, 1.0);
    return TapeVariable(tape, idx, val);
}

inline TapeVariable operator-(const TapeVariable& a, const TapeVariable& b) {
    Tape* tape = a.tape_ ? a.tape_ : b.tape_;
    if (!tape)
        return TapeVariable(nullptr, -1, a.value_ - b.value_);

    double val = a.value_ - b.value_;
    int idx =
        tape->RecordOperation(OpCode::kSub, val, a.index_, b.index_, 1.0, -1.0);
    return TapeVariable(tape, idx, val);
}

inline TapeVariable operator*(const TapeVariable& a, const TapeVariable& b) {
    Tape* tape = a.tape_ ? a.tape_ : b.tape_;
    if (!tape)
        return TapeVariable(nullptr, -1, a.value_ * b.value_);

    double val = a.value_ * b.value_;
    // d(a*b)/da = b, d(a*b)/db = a
    int idx = tape->RecordOperation(OpCode::kMul, val, a.index_, b.index_,
                                    b.value_, a.value_);
    return TapeVariable(tape, idx, val);
}

inline TapeVariable operator/(const TapeVariable& a, const TapeVariable& b) {
    Tape* tape = a.tape_ ? a.tape_ : b.tape_;
    if (!tape)
        return TapeVariable(nullptr, -1, a.value_ / b.value_);

    double val = a.value_ / b.value_;
    // d(a/b)/da = 1/b, d(a/b)/db = -a/b²
    double inv_b = 1.0 / b.value_;
    int idx = tape->RecordOperation(OpCode::kDiv, val, a.index_, b.index_,
                                    inv_b, -a.value_ * inv_b * inv_b);
    return TapeVariable(tape, idx, val);
}

// Scalar-variable operators

inline TapeVariable operator+(double s, const TapeVariable& v) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, s + v.value_);
    TapeVariable sc = v.tape_->CreateConstant(s);
    return sc + v;
}

inline TapeVariable operator-(double s, const TapeVariable& v) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, s - v.value_);
    TapeVariable sc = v.tape_->CreateConstant(s);
    return sc - v;
}

inline TapeVariable operator*(double s, const TapeVariable& v) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, s * v.value_);
    double val = s * v.value_;
    int idx = v.tape_->RecordOperation(OpCode::kMul, val, v.index_, -1, s, 0.0);
    return TapeVariable(v.tape_, idx, val);
}

inline TapeVariable operator/(double s, const TapeVariable& v) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, s / v.value_);
    TapeVariable sc = v.tape_->CreateConstant(s);
    return sc / v;
}

inline TapeVariable operator+(const TapeVariable& v, double s) {
    return s + v;
}

inline TapeVariable operator-(const TapeVariable& v, double s) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, v.value_ - s);
    TapeVariable sc = v.tape_->CreateConstant(s);
    return v - sc;
}

inline TapeVariable operator*(const TapeVariable& v, double s) {
    return s * v;
}

inline TapeVariable operator/(const TapeVariable& v, double s) {
    if (!v.tape_)
        return TapeVariable(nullptr, -1, v.value_ / s);
    double val = v.value_ / s;
    int idx =
        v.tape_->RecordOperation(OpCode::kDiv, val, v.index_, -1, 1.0 / s, 0.0);
    return TapeVariable(v.tape_, idx, val);
}

// Mathematical functions

inline TapeVariable sqrt(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::sqrt(x.value_));
    double val = std::sqrt(x.value_);
    // d/dx sqrt(x) = 0.5/sqrt(x)
    int idx = x.tape_->RecordOperation(OpCode::kSqrt, val, x.index_, -1,
                                       0.5 / val, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable sin(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::sin(x.value_));
    double val = std::sin(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kSin, val, x.index_, -1,
                                       std::cos(x.value_), 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable cos(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::cos(x.value_));
    double val = std::cos(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kCos, val, x.index_, -1,
                                       -std::sin(x.value_), 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable tan(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::tan(x.value_));
    double val = std::tan(x.value_);
    double c = std::cos(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kTan, val, x.index_, -1,
                                       1.0 / (c * c), 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable exp(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::exp(x.value_));
    double val = std::exp(x.value_);
    int idx =
        x.tape_->RecordOperation(OpCode::kExp, val, x.index_, -1, val, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable log(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::log(x.value_));
    double val = std::log(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kLog, val, x.index_, -1,
                                       1.0 / x.value_, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable pow(const TapeVariable& base,
                        const TapeVariable& exponent) {
    Tape* tape = base.tape_ ? base.tape_ : exponent.tape_;
    if (!tape) {
        return TapeVariable(nullptr, -1,
                            std::pow(base.value_, exponent.value_));
    }
    double val = std::pow(base.value_, exponent.value_);
    // d/da a^b = b * a^(b-1), d/db a^b = a^b * log(a)
    double p1 = exponent.value_ * std::pow(base.value_, exponent.value_ - 1.0);
    double p2 = val * std::log(base.value_);
    int idx = tape->RecordOperation(OpCode::kPow, val, base.index_,
                                    exponent.index_, p1, p2);
    return TapeVariable(tape, idx, val);
}

inline TapeVariable pow(const TapeVariable& base, double exponent) {
    if (!base.tape_) {
        return TapeVariable(nullptr, -1, std::pow(base.value_, exponent));
    }
    double val = std::pow(base.value_, exponent);
    double p1 = exponent * std::pow(base.value_, exponent - 1.0);
    int idx = base.tape_->RecordOperation(OpCode::kPow, val, base.index_, -1,
                                          p1, 0.0);
    return TapeVariable(base.tape_, idx, val);
}

inline TapeVariable pow(double base, const TapeVariable& exponent) {
    if (!exponent.tape_) {
        return TapeVariable(nullptr, -1, std::pow(base, exponent.value_));
    }
    double val = std::pow(base, exponent.value_);
    double p2 = val * std::log(base);
    int idx = exponent.tape_->RecordOperation(OpCode::kPow, val, -1,
                                              exponent.index_, 0.0, p2);
    return TapeVariable(exponent.tape_, idx, val);
}

inline TapeVariable abs(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::abs(x.value_));
    double val = std::abs(x.value_);
    double sign = (x.value_ >= 0) ? 1.0 : -1.0;
    int idx =
        x.tape_->RecordOperation(OpCode::kAbs, val, x.index_, -1, sign, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable asin(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::asin(x.value_));
    double val = std::asin(x.value_);
    double p1 = 1.0 / std::sqrt(1.0 - x.value_ * x.value_);
    int idx =
        x.tape_->RecordOperation(OpCode::kAsin, val, x.index_, -1, p1, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable acos(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::acos(x.value_));
    double val = std::acos(x.value_);
    double p1 = -1.0 / std::sqrt(1.0 - x.value_ * x.value_);
    int idx =
        x.tape_->RecordOperation(OpCode::kAcos, val, x.index_, -1, p1, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable atan(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::atan(x.value_));
    double val = std::atan(x.value_);
    double p1 = 1.0 / (1.0 + x.value_ * x.value_);
    int idx =
        x.tape_->RecordOperation(OpCode::kAtan, val, x.index_, -1, p1, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable atan2(const TapeVariable& y, const TapeVariable& x) {
    Tape* tape = y.tape_ ? y.tape_ : x.tape_;
    if (!tape) {
        return TapeVariable(nullptr, -1, std::atan2(y.value_, x.value_));
    }
    double val = std::atan2(y.value_, x.value_);
    double r2 = x.value_ * x.value_ + y.value_ * y.value_;
    double p1 = x.value_ / r2;   // d/dy atan2(y,x)
    double p2 = -y.value_ / r2;  // d/dx atan2(y,x)
    int idx =
        tape->RecordOperation(OpCode::kAtan2, val, y.index_, x.index_, p1, p2);
    return TapeVariable(tape, idx, val);
}

inline TapeVariable sinh(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::sinh(x.value_));
    double val = std::sinh(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kSinh, val, x.index_, -1,
                                       std::cosh(x.value_), 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable cosh(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::cosh(x.value_));
    double val = std::cosh(x.value_);
    int idx = x.tape_->RecordOperation(OpCode::kCosh, val, x.index_, -1,
                                       std::sinh(x.value_), 0.0);
    return TapeVariable(x.tape_, idx, val);
}

inline TapeVariable tanh(const TapeVariable& x) {
    if (!x.tape_)
        return TapeVariable(nullptr, -1, std::tanh(x.value_));
    double val = std::tanh(x.value_);
    double sech2 = 1.0 - val * val;
    int idx =
        x.tape_->RecordOperation(OpCode::kTanh, val, x.index_, -1, sech2, 0.0);
    return TapeVariable(x.tape_, idx, val);
}

// Convenience functions for using reverse-mode AD

/// Compute gradient of a scalar function using reverse-mode AD
///
/// @tparam F Function type: takes std::vector<TapeVariable>& and returns
/// TapeVariable
/// @param f The function to differentiate
/// @param x The point at which to compute the gradient
/// @return Gradient vector at x
template <typename F>
std::vector<double> ComputeGradient(F&& f, const std::vector<double>& x) {
    Tape tape;
    tape.Activate();

    std::vector<TapeVariable> vars;
    vars.reserve(x.size());
    for (double xi : x) {
        vars.push_back(tape.CreateVariable(xi));
    }

    TapeVariable result = f(vars);

    auto gradient = tape.ComputeGradient(result);

    tape.Deactivate();
    return gradient;
}

/// RAII helper to activate/deactivate a tape
class TapeScope {
public:
    explicit TapeScope(Tape& tape) : tape_(tape) { tape_.Activate(); }
    ~TapeScope() { tape_.Deactivate(); }

    TapeScope(const TapeScope&) = delete;
    TapeScope& operator=(const TapeScope&) = delete;

private:
    Tape& tape_;
};

}  // namespace gradients
}  // namespace pegasus

#endif  // PEGASUS_GRADIENTS_AUTODIFF_TAPE_H_
