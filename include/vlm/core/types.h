// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_CORE_TYPES_H_
#define VLM_CORE_TYPES_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cstdint>

namespace vlm {

/// @defgroup CoreTypes Core Type Definitions
/// @brief Fundamental types and aliases used throughout the VLM library.
/// @{

// Scalar type - can be changed to float for single precision
using Scalar = double;

// Fixed-size vector types
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;

// Fixed-size matrix types
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;
using Mat4 = Eigen::Matrix4d;

// Dynamic-size types
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;

// Row vector types (for compatibility with matrix operations)
using RowVec3 = Eigen::RowVector3d;
using RowVecX = Eigen::RowVectorXd;

// Index type for Eigen compatibility
using Index = Eigen::Index;

/// @}

/// @defgroup Constants Physical and Mathematical Constants
/// @{

/// Mathematical constant pi
constexpr Scalar kPi = 3.14159265358979323846;

/// Conversion factor from degrees to radians
constexpr Scalar kDegreesToRadians = kPi / 180.0;

/// Conversion factor from radians to degrees
constexpr Scalar kRadiansToDegrees = 180.0 / kPi;

/// Small value for numerical comparisons
constexpr Scalar kEpsilon = 1.0e-12;

/// Default tolerance for iterative solvers
constexpr Scalar kDefaultTolerance = 1.0e-10;

/// @}

/// @brief Symmetry type for XZ and XY planes.
///
/// Controls how the geometry is mirrored across symmetry planes.
enum class SymmetryType : int {
  kNone = 0,           ///< No symmetry
  kSymmetric = 1,      ///< Regular symmetry (same sign)
  kAntisymmetric = -1  ///< Antisymmetric (opposite sign, free surface)
};

/// @brief Spacing distribution type for panel generation.
enum class SpacingType {
  kUniform,   ///< Uniform spacing
  kCosine,    ///< Cosine spacing (clustered at edges)
  kSine,      ///< Sine spacing
  kBlend      ///< Blend of uniform and cosine
};

}  // namespace vlm

#endif  // VLM_CORE_TYPES_H_
