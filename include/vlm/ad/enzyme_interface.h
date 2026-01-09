// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file enzyme_interface.h
/// @brief Enzyme automatic differentiation interface declarations
///
/// This header provides the interface to Enzyme AD when compiled with
/// VLM_USE_ENZYME. Enzyme operates at the LLVM IR level, automatically
/// differentiating compiled code without requiring source modifications.
///
/// To enable Enzyme, configure with:
///   cmake -DVLM_ENABLE_ENZYME=ON -DCMAKE_CXX_COMPILER=clang++ ..
///
/// Enzyme must be installed (e.g., `brew install enzyme` on macOS).

#ifdef VLM_USE_ENZYME

namespace vlm {
namespace ad {

// ============================================================================
// Enzyme External Declarations
// ============================================================================

/// Activity annotations for Enzyme
/// These are used to mark how arguments should be treated during differentiation
extern "C" {
/// Marker for duplicated arguments (active variables with shadow)
extern int enzyme_dup;
/// Marker for constant arguments (not differentiated)
extern int enzyme_const;
/// Marker for output-only arguments
extern int enzyme_out;
/// Marker for duplicated pointer arguments without need for primal
extern int enzyme_dupnoneed;
}

// ============================================================================
// Enzyme Autodiff Function Templates
// ============================================================================

/// @brief Reverse-mode automatic differentiation (adjoint)
///
/// Computes gradients of a scalar output with respect to all active inputs.
/// This is efficient when there are many inputs and few outputs (e.g., optimization).
///
/// Usage:
/// @code
/// double f(double x, double y) { return x*x + y*y; }
///
/// double x = 3.0, y = 4.0;
/// double dx = 0.0, dy = 0.0;
/// __enzyme_autodiff<void>((void*)f,
///     enzyme_dup, x, &dx,    // x is active, gradient stored in dx
///     enzyme_dup, y, &dy);   // y is active, gradient stored in dy
/// // Now dx = 6.0 (∂f/∂x), dy = 8.0 (∂f/∂y)
/// @endcode
///
/// @tparam Return Return type of the original function
/// @tparam Args Argument types (alternating: annotation, value, shadow pairs)
/// @param fn Pointer to the function to differentiate
/// @param args Activity annotations and argument/shadow pairs
/// @return For scalar functions, returns the gradient; otherwise void
template <typename Return, typename... Args>
Return __enzyme_autodiff(void* fn, Args... args);

/// @brief Forward-mode automatic differentiation (tangent)
///
/// Computes directional derivatives of all outputs with respect to a single
/// input perturbation. Efficient when there are few inputs and many outputs.
///
/// Usage:
/// @code
/// void f(double x, double* out) { out[0] = x*x; out[1] = x*x*x; }
///
/// double x = 2.0, dx = 1.0;  // Seed dx=1 to get d/dx
/// double out[2], dout[2];
/// __enzyme_fwddiff<void>((void*)f,
///     enzyme_dup, x, dx,           // x with tangent dx
///     enzyme_dup, out, dout);      // outputs with tangent outputs
/// // Now dout[0] = 4.0 (d(x²)/dx), dout[1] = 12.0 (d(x³)/dx)
/// @endcode
///
/// @tparam Return Return type of the original function
/// @tparam Args Argument types (alternating: annotation, value, tangent pairs)
/// @param fn Pointer to the function to differentiate
/// @param args Activity annotations and argument/tangent pairs
/// @return The function return value (if non-void) plus tangent
template <typename Return, typename... Args>
Return __enzyme_fwddiff(void* fn, Args... args);

// ============================================================================
// Utility Functions
// ============================================================================

/// @brief Check if Enzyme AD is available at runtime
/// @return true if compiled with Enzyme support
constexpr bool IsEnzymeAvailable() { return true; }

// ============================================================================
// Helper Macros
// ============================================================================

/// Helper macro to mark a variable as active (duplicated) in autodiff calls
#define VLM_ENZYME_DUP enzyme_dup

/// Helper macro to mark a variable as constant in autodiff calls
#define VLM_ENZYME_CONST enzyme_const

/// Helper macro to mark a variable as output-only in autodiff calls
#define VLM_ENZYME_OUT enzyme_out

}  // namespace ad
}  // namespace vlm

#else  // !VLM_USE_ENZYME

namespace vlm {
namespace ad {

/// @brief Check if Enzyme AD is available at runtime
/// @return false when compiled without Enzyme support
constexpr bool IsEnzymeAvailable() { return false; }

// Stub macros when Enzyme is not available
#define VLM_ENZYME_DUP 0
#define VLM_ENZYME_CONST 0
#define VLM_ENZYME_OUT 0

}  // namespace ad
}  // namespace vlm

#endif  // VLM_USE_ENZYME
