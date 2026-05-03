/*
    Project VLM
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

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

#ifndef VLM_UTILS_MATH_UTILITIES_H_
#define VLM_UTILS_MATH_UTILITIES_H_

#include "vlm/core/types.h"

#include <cmath>

namespace vlm {

/// @defgroup MathUtilities Mathematical Utility Functions
/// @brief Helper functions for aerodynamic calculations.
/// @{

/// @brief Computes the Prandtl-Glauert compressibility factor.
///
/// The Prandtl-Glauert factor beta = sqrt(1 - M^2) is used to correct
/// for compressibility effects in subsonic flow.
///
/// @param mach_number Freestream Mach number (must be < 1)
/// @return Prandtl-Glauert factor beta
inline Scalar ComputePrandtlGlauertFactor(Scalar mach_number) {
    return std::sqrt(1.0 - mach_number * mach_number);
}

/// @brief Computes induced velocity from a finite vortex segment.
///
/// Uses the Biot-Savart law with finite-core regularization to compute
/// the velocity induced at a field point by a straight vortex segment
/// of unit circulation.
///
/// The velocity is computed as:
///   V = (Gamma / 4*pi) * [(r1 x r2) / |r1 x r2|^2] * (r0 . (r1/|r1| -
///   r2/|r2|))
///
/// where r1 and r2 are vectors from segment endpoints to the field point,
/// and r0 is the segment vector.
///
/// The finite core prevents singularities when the field point approaches
/// the vortex line. The core radius should be a small fraction of the
/// local panel dimension.
///
/// @param field_point Point where velocity is computed
/// @param segment_start Start point of vortex segment
/// @param segment_end End point of vortex segment
/// @param core_radius Finite core radius for regularization
/// @return Induced velocity vector (per unit circulation)
///
/// Reference: Katz & Plotkin, "Low-Speed Aerodynamics", Section 2.5
inline Vec3 ComputeVortexSegmentVelocity(const Vec3& field_point,
                                         const Vec3& segment_start,
                                         const Vec3& segment_end,
                                         Scalar core_radius) {
    // Vectors from segment endpoints to field point
    Vec3 r1 = field_point - segment_start;
    Vec3 r2 = field_point - segment_end;

    // Segment vector
    Vec3 r0 = segment_end - segment_start;

    // Cross product r1 x r2
    Vec3 cross = r1.cross(r2);
    Scalar cross_mag_sq = cross.squaredNorm();

    // Regularization: add core radius squared to prevent singularity
    Scalar core_sq = core_radius * core_radius;
    cross_mag_sq = std::max(cross_mag_sq, core_sq);

    // Magnitudes
    Scalar r1_mag = r1.norm();
    Scalar r2_mag = r2.norm();

    // Avoid division by zero for degenerate cases
    if (r1_mag < kEpsilon || r2_mag < kEpsilon) {
        return Vec3::Zero();
    }

    // Biot-Savart formula
    Scalar factor =
        r0.dot(r1 / r1_mag - r2 / r2_mag) / (4.0 * kPi * cross_mag_sq);

    return factor * cross;
}

/// @brief Computes induced velocity from a semi-infinite trailing vortex.
///
/// Computes the velocity induced by a semi-infinite vortex starting at
/// a given point and extending in the downstream direction (positive x).
///
/// @param field_point Point where velocity is computed
/// @param vortex_start Starting point of the semi-infinite vortex
/// @param core_radius Finite core radius for regularization
/// @return Induced velocity vector (per unit circulation)
inline Vec3 ComputeTrailingVortexVelocity(const Vec3& field_point,
                                          const Vec3& vortex_start,
                                          Scalar core_radius) {
    // Vector from vortex start to field point
    Vec3 r = field_point - vortex_start;

    // For a semi-infinite vortex along +x, the cross product with unit x is:
    // (1,0,0) x r = (0, r.z, -r.y)
    Vec3 cross(0.0, r.z(), -r.y());
    Scalar cross_mag_sq = cross.squaredNorm();

    // Regularization
    Scalar core_sq = core_radius * core_radius;
    cross_mag_sq = std::max(cross_mag_sq, core_sq);

    // Distance along x-axis
    Scalar x = r.x();
    Scalar r_mag = r.norm();

    if (r_mag < kEpsilon) {
        return Vec3::Zero();
    }

    // Biot-Savart for semi-infinite vortex
    // V = Gamma / (4*pi) * (1 + x/|r|) * (x_hat x r) / |x_hat x r|^2
    Scalar factor = (1.0 + x / r_mag) / (4.0 * kPi * cross_mag_sq);

    return factor * cross;
}

/// @brief Generates cosine-spaced points in the range [0, 1].
///
/// Cosine spacing concentrates points near the endpoints, which is useful
/// for capturing leading and trailing edge effects in aerodynamic panels.
///
/// @param num_points Number of points to generate (must be >= 2)
/// @param values Output vector of cosine-spaced values
inline void GenerateCosineSpacing(Index num_points,
                                  std::vector<Scalar>& values) {
    values.resize(num_points);
    for (Index i = 0; i < num_points; ++i) {
        Scalar theta =
            kPi * static_cast<Scalar>(i) / static_cast<Scalar>(num_points - 1);
        values[i] = 0.5 * (1.0 - std::cos(theta));
    }
}

/// @brief Generates uniform-spaced points in the range [0, 1].
///
/// @param num_points Number of points to generate (must be >= 2)
/// @param values Output vector of uniformly-spaced values
inline void GenerateUniformSpacing(Index num_points,
                                   std::vector<Scalar>& values) {
    values.resize(num_points);
    for (Index i = 0; i < num_points; ++i) {
        values[i] =
            static_cast<Scalar>(i) / static_cast<Scalar>(num_points - 1);
    }
}

/// @brief Generates sine-spaced points in the range [0, 1].
///
/// Sine spacing concentrates points near the trailing edge (x=1).
///
/// @param num_points Number of points to generate (must be >= 2)
/// @param values Output vector of sine-spaced values
inline void GenerateSineSpacing(Index num_points, std::vector<Scalar>& values) {
    values.resize(num_points);
    for (Index i = 0; i < num_points; ++i) {
        Scalar theta = 0.5 * kPi * static_cast<Scalar>(i) /
                       static_cast<Scalar>(num_points - 1);
        values[i] = std::sin(theta);
    }
}

/// @brief Blends between uniform and cosine spacing.
///
/// @param num_points Number of points to generate
/// @param blend_factor Blend factor: 0 = uniform, 1 = full cosine
/// @param values Output vector of blended spacing values
inline void GenerateBlendedSpacing(Index num_points, Scalar blend_factor,
                                   std::vector<Scalar>& values) {
    std::vector<Scalar> uniform, cosine;
    GenerateUniformSpacing(num_points, uniform);
    GenerateCosineSpacing(num_points, cosine);

    values.resize(num_points);
    for (Index i = 0; i < num_points; ++i) {
        values[i] =
            (1.0 - blend_factor) * uniform[i] + blend_factor * cosine[i];
    }
}

/// @brief Generates spacing points based on spacing type.
///
/// @param num_points Number of points to generate (must be >= 2)
/// @param spacing_type Type of spacing distribution
/// @return Vector of spacing values in range [0, 1]
inline std::vector<Scalar> GenerateCosineSpacing(Index num_points,
                                                 SpacingType spacing_type) {
    std::vector<Scalar> values;
    switch (spacing_type) {
        case SpacingType::kUniform:
        case SpacingType::kNone:
            GenerateUniformSpacing(num_points, values);
            break;
        case SpacingType::kCosine:
            GenerateCosineSpacing(num_points, values);
            break;
        case SpacingType::kSine:
            GenerateSineSpacing(num_points, values);
            break;
        case SpacingType::kBlend:
            GenerateBlendedSpacing(num_points, 0.5, values);
            break;
    }
    return values;
}

/// @brief Rotates a vector about the x-axis.
///
/// @param v Vector to rotate
/// @param angle Rotation angle in radians
/// @return Rotated vector
inline Vec3 RotateAboutX(const Vec3& v, Scalar angle) {
    Scalar c = std::cos(angle);
    Scalar s = std::sin(angle);
    return {v.x(), c * v.y() - s * v.z(), s * v.y() + c * v.z()};
}

/// @brief Rotates a vector about the y-axis.
///
/// @param v Vector to rotate
/// @param angle Rotation angle in radians
/// @return Rotated vector
inline Vec3 RotateAboutY(const Vec3& v, Scalar angle) {
    Scalar c = std::cos(angle);
    Scalar s = std::sin(angle);
    return {c * v.x() + s * v.z(), v.y(), -s * v.x() + c * v.z()};
}

/// @brief Rotates a vector about the z-axis.
///
/// @param v Vector to rotate
/// @param angle Rotation angle in radians
/// @return Rotated vector
inline Vec3 RotateAboutZ(const Vec3& v, Scalar angle) {
    Scalar c = std::cos(angle);
    Scalar s = std::sin(angle);
    return {c * v.x() - s * v.y(), s * v.x() + c * v.y(), v.z()};
}

/// @brief Computes rotation matrix for wind-to-body axis transformation.
///
/// Transforms from wind axes (aligned with freestream) to body axes
/// using angle of attack (alpha) and sideslip angle (beta).
///
/// @param alpha Angle of attack in radians
/// @param beta Sideslip angle in radians
/// @return 3x3 rotation matrix from wind to body axes
inline Mat3 WindToBodyRotation(Scalar alpha, Scalar beta) {
    Scalar ca = std::cos(alpha);
    Scalar sa = std::sin(alpha);
    Scalar cb = std::cos(beta);
    Scalar sb = std::sin(beta);

    Mat3 rotation;
    rotation << ca * cb, -ca * sb, -sa, sb, cb, 0.0, sa * cb, -sa * sb, ca;
    return rotation;
}

/// @brief Computes rotation matrix for body-to-stability axis transformation.
///
/// Stability axes are rotated from body axes by the angle of attack only.
///
/// @param alpha Angle of attack in radians
/// @return 3x3 rotation matrix from body to stability axes
inline Mat3 BodyToStabilityRotation(Scalar alpha) {
    Scalar ca = std::cos(alpha);
    Scalar sa = std::sin(alpha);

    Mat3 rotation;
    rotation << ca, 0.0, sa, 0.0, 1.0, 0.0, -sa, 0.0, ca;
    return rotation;
}

/// @}

}  // namespace vlm

#endif  // VLM_UTILS_MATH_UTILITIES_H_
