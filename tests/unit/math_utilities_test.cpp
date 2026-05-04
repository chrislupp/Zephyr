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

#include "pegasus/utils/math_utilities.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTol = 1.0e-10;

// ============================================================================
// Prandtl-Glauert Factor
// ============================================================================

TEST(MathUtilitiesTest, PrandtlGlauertZeroMach) {
    EXPECT_NEAR(ComputePrandtlGlauertFactor(0.0), 1.0, kTol);
}

TEST(MathUtilitiesTest, PrandtlGlauertSubsonic) {
    double beta = ComputePrandtlGlauertFactor(0.5);
    EXPECT_NEAR(beta, std::sqrt(1.0 - 0.25), kTol);
}

// ============================================================================
// Vortex Segment Velocity
// ============================================================================

TEST(MathUtilitiesTest, VortexSegmentAboveCenter) {
    Vec3 field(0.5, 0.0, 1.0);
    Vec3 start(0.0, 0.0, 0.0);
    Vec3 end(1.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeVortexSegmentVelocity(field, start, end, core);
    EXPECT_TRUE(vel.allFinite());
    // Velocity should be in the y-direction (by right-hand rule)
    EXPECT_GT(std::abs(vel.y()), 1e-6);
}

TEST(MathUtilitiesTest, VortexSegmentOnAxis) {
    // Field point on the vortex axis - should return zero (regularized)
    Vec3 field(0.5, 0.0, 0.0);
    Vec3 start(0.0, 0.0, 0.0);
    Vec3 end(1.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeVortexSegmentVelocity(field, start, end, core);
    EXPECT_TRUE(vel.allFinite());
}

TEST(MathUtilitiesTest, VortexSegmentAtEndpoint) {
    // Field point at segment endpoint - degenerate case
    Vec3 field(0.0, 0.0, 0.0);
    Vec3 start(0.0, 0.0, 0.0);
    Vec3 end(1.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeVortexSegmentVelocity(field, start, end, core);
    EXPECT_TRUE(vel.allFinite());
    EXPECT_NEAR(vel.norm(), 0.0, 1e-6);
}

// ============================================================================
// Trailing Vortex Velocity
// ============================================================================

TEST(MathUtilitiesTest, TrailingVortexAbove) {
    Vec3 field(0.0, 0.0, 1.0);
    Vec3 start(0.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeTrailingVortexVelocity(field, start, core);
    EXPECT_TRUE(vel.allFinite());
    // cross = (0, r.z, -r.y) = (0, 1, 0) -> positive y velocity
    EXPECT_GT(vel.y(), 0.0);
}

TEST(MathUtilitiesTest, TrailingVortexAtStart) {
    Vec3 field(0.0, 0.0, 0.0);
    Vec3 start(0.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeTrailingVortexVelocity(field, start, core);
    EXPECT_TRUE(vel.allFinite());
    EXPECT_NEAR(vel.norm(), 0.0, 1e-6);
}

TEST(MathUtilitiesTest, TrailingVortexDownstream) {
    Vec3 field(5.0, 0.0, 1.0);
    Vec3 start(0.0, 0.0, 0.0);
    double core = 0.01;

    Vec3 vel = ComputeTrailingVortexVelocity(field, start, core);
    EXPECT_TRUE(vel.allFinite());
    EXPECT_GT(vel.y(), 0.0);
}

// ============================================================================
// Spacing Functions
// ============================================================================

TEST(MathUtilitiesTest, CosineSpacingEndpoints) {
    std::vector<Scalar> values;
    GenerateCosineSpacing(11, values);

    EXPECT_EQ(values.size(), 11u);
    EXPECT_NEAR(values.front(), 0.0, kTol);
    EXPECT_NEAR(values.back(), 1.0, kTol);
}

TEST(MathUtilitiesTest, CosineSpacingConcentration) {
    std::vector<Scalar> values;
    GenerateCosineSpacing(21, values);

    // Points should be denser near 0 and 1
    double first_gap = values[1] - values[0];
    double mid_gap = values[11] - values[10];
    EXPECT_LT(first_gap, mid_gap);
}

TEST(MathUtilitiesTest, UniformSpacing) {
    std::vector<Scalar> values;
    GenerateUniformSpacing(5, values);

    EXPECT_EQ(values.size(), 5u);
    EXPECT_NEAR(values[0], 0.0, kTol);
    EXPECT_NEAR(values[1], 0.25, kTol);
    EXPECT_NEAR(values[2], 0.5, kTol);
    EXPECT_NEAR(values[3], 0.75, kTol);
    EXPECT_NEAR(values[4], 1.0, kTol);
}

TEST(MathUtilitiesTest, SineSpacing) {
    std::vector<Scalar> values;
    GenerateSineSpacing(11, values);

    EXPECT_EQ(values.size(), 11u);
    EXPECT_NEAR(values.front(), 0.0, kTol);
    EXPECT_NEAR(values.back(), 1.0, kTol);

    // Sine spacing: denser near 0
    double first_gap = values[1] - values[0];
    double last_gap = values[10] - values[9];
    EXPECT_GT(first_gap, last_gap);
}

TEST(MathUtilitiesTest, BlendedSpacingUniform) {
    // blend_factor = 0 should give uniform
    std::vector<Scalar> blended, uniform;
    GenerateBlendedSpacing(5, 0.0, blended);
    GenerateUniformSpacing(5, uniform);

    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR(blended[i], uniform[i], kTol);
    }
}

TEST(MathUtilitiesTest, BlendedSpacingCosine) {
    // blend_factor = 1 should give cosine
    std::vector<Scalar> blended, cosine;
    GenerateBlendedSpacing(5, 1.0, blended);
    GenerateCosineSpacing(5, cosine);

    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR(blended[i], cosine[i], kTol);
    }
}

TEST(MathUtilitiesTest, SpacingDispatcherUniform) {
    auto values = GenerateCosineSpacing(5, SpacingType::kUniform);
    EXPECT_EQ(values.size(), 5u);
    EXPECT_NEAR(values[2], 0.5, kTol);
}

TEST(MathUtilitiesTest, SpacingDispatcherCosine) {
    auto values = GenerateCosineSpacing(5, SpacingType::kCosine);
    EXPECT_EQ(values.size(), 5u);
    EXPECT_NEAR(values.front(), 0.0, kTol);
    EXPECT_NEAR(values.back(), 1.0, kTol);
}

TEST(MathUtilitiesTest, SpacingDispatcherSine) {
    auto values = GenerateCosineSpacing(5, SpacingType::kSine);
    EXPECT_EQ(values.size(), 5u);
}

TEST(MathUtilitiesTest, SpacingDispatcherBlend) {
    auto values = GenerateCosineSpacing(5, SpacingType::kBlend);
    EXPECT_EQ(values.size(), 5u);
}

TEST(MathUtilitiesTest, SpacingDispatcherNone) {
    auto values = GenerateCosineSpacing(5, SpacingType::kNone);
    EXPECT_EQ(values.size(), 5u);
    // kNone should give uniform
    EXPECT_NEAR(values[2], 0.5, kTol);
}

// ============================================================================
// Rotation Functions
// ============================================================================

TEST(MathUtilitiesTest, RotateAboutXIdentity) {
    Vec3 v(1.0, 2.0, 3.0);
    Vec3 r = RotateAboutX(v, 0.0);
    EXPECT_NEAR(r.x(), v.x(), kTol);
    EXPECT_NEAR(r.y(), v.y(), kTol);
    EXPECT_NEAR(r.z(), v.z(), kTol);
}

TEST(MathUtilitiesTest, RotateAboutX90) {
    Vec3 v(0.0, 1.0, 0.0);
    Vec3 r = RotateAboutX(v, kPi / 2.0);
    EXPECT_NEAR(r.x(), 0.0, kTol);
    EXPECT_NEAR(r.y(), 0.0, kTol);
    EXPECT_NEAR(r.z(), 1.0, kTol);
}

TEST(MathUtilitiesTest, RotateAboutYIdentity) {
    Vec3 v(1.0, 2.0, 3.0);
    Vec3 r = RotateAboutY(v, 0.0);
    EXPECT_NEAR(r.x(), v.x(), kTol);
    EXPECT_NEAR(r.y(), v.y(), kTol);
    EXPECT_NEAR(r.z(), v.z(), kTol);
}

TEST(MathUtilitiesTest, RotateAboutY90) {
    Vec3 v(1.0, 0.0, 0.0);
    Vec3 r = RotateAboutY(v, kPi / 2.0);
    EXPECT_NEAR(r.x(), 0.0, kTol);
    EXPECT_NEAR(r.y(), 0.0, kTol);
    EXPECT_NEAR(r.z(), -1.0, kTol);
}

TEST(MathUtilitiesTest, RotateAboutZIdentity) {
    Vec3 v(1.0, 2.0, 3.0);
    Vec3 r = RotateAboutZ(v, 0.0);
    EXPECT_NEAR(r.x(), v.x(), kTol);
    EXPECT_NEAR(r.y(), v.y(), kTol);
    EXPECT_NEAR(r.z(), v.z(), kTol);
}

TEST(MathUtilitiesTest, RotateAboutZ90) {
    Vec3 v(1.0, 0.0, 0.0);
    Vec3 r = RotateAboutZ(v, kPi / 2.0);
    EXPECT_NEAR(r.x(), 0.0, kTol);
    EXPECT_NEAR(r.y(), 1.0, kTol);
    EXPECT_NEAR(r.z(), 0.0, kTol);
}

// ============================================================================
// Coordinate Transformations
// ============================================================================

TEST(MathUtilitiesTest, WindToBodyIdentity) {
    Mat3 R = WindToBodyRotation(0.0, 0.0);
    Mat3 I = Mat3::Identity();
    EXPECT_TRUE((R - I).norm() < kTol);
}

TEST(MathUtilitiesTest, WindToBodyOrthogonality) {
    Mat3 R = WindToBodyRotation(0.1, 0.05);
    Mat3 RtR = R.transpose() * R;
    Mat3 I = Mat3::Identity();
    EXPECT_LT((RtR - I).norm(), 1e-10);
}

TEST(MathUtilitiesTest, WindToBodyAlphaOnly) {
    double alpha = kPi / 6.0;
    Mat3 R = WindToBodyRotation(alpha, 0.0);
    // First row: [cos(a), 0, -sin(a)]
    EXPECT_NEAR(R(0, 0), std::cos(alpha), kTol);
    EXPECT_NEAR(R(0, 2), -std::sin(alpha), kTol);
}

TEST(MathUtilitiesTest, BodyToStabilityIdentity) {
    Mat3 R = BodyToStabilityRotation(0.0);
    Mat3 I = Mat3::Identity();
    EXPECT_TRUE((R - I).norm() < kTol);
}

TEST(MathUtilitiesTest, BodyToStabilityOrthogonality) {
    Mat3 R = BodyToStabilityRotation(0.2);
    Mat3 RtR = R.transpose() * R;
    Mat3 I = Mat3::Identity();
    EXPECT_LT((RtR - I).norm(), 1e-10);
}

TEST(MathUtilitiesTest, BodyToStabilityElements) {
    double alpha = kPi / 4.0;
    Mat3 R = BodyToStabilityRotation(alpha);
    EXPECT_NEAR(R(0, 0), std::cos(alpha), kTol);
    EXPECT_NEAR(R(0, 2), std::sin(alpha), kTol);
    EXPECT_NEAR(R(1, 1), 1.0, kTol);
    EXPECT_NEAR(R(2, 0), -std::sin(alpha), kTol);
    EXPECT_NEAR(R(2, 2), std::cos(alpha), kTol);
}

}  // namespace
}  // namespace pegasus
