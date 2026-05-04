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

#include "pegasus/gradients/analytical_trefftz.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace pegasus::gradients {
namespace {

constexpr double kFDStep = 1e-7;
constexpr double kFDTol = 1e-4;

class TrefftzGradientTest : public ::testing::Test {
protected:
    void SetUp() override {
        num_strips_ = 3;
        ref_area_ = 10.0;
        ref_span_ = 10.0;
        core_chord_frac_ = 0.25;
        core_width_frac_ = 0.05;

        // 3-strip wing: left TE, right TE, with some y-span
        gamma_ = {0.10, 0.12, 0.08};
        // TE left points: x doesn't matter for Trefftz, only y,z
        te_left_ = {1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 1.0, 1.0, 0.0};
        te_right_ = {1.0, 0.5, 0.0, 1.0, 1.0, 0.0, 1.0, 1.5, 0.0};
        chord_ = {1.0, 1.0, 1.0};
        surface_idx_ = {0, 0, 0};
    }

    int num_strips_;
    double ref_area_, ref_span_;
    double core_chord_frac_, core_width_frac_;
    std::vector<double> gamma_, te_left_, te_right_, chord_;
    std::vector<int> surface_idx_;
};

TEST_F(TrefftzGradientTest, PrimalKernelComputes) {
    double CL, CY, CDi;
    compute_trefftz_drag_kernel(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, &CL, &CY, &CDi);

    EXPECT_GT(CL, 0.0);
    EXPECT_NEAR(CY, 0.0, 1e-10);
    EXPECT_GT(CDi, 0.0);
    EXPECT_TRUE(std::isfinite(CL));
    EXPECT_TRUE(std::isfinite(CDi));
}

TEST_F(TrefftzGradientTest, SpanEfficiency) {
    double CL, CY, CDi;
    compute_trefftz_drag_kernel(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, &CL, &CY, &CDi);

    double e = compute_span_efficiency(CL, CY, CDi, ref_area_, ref_span_);
    EXPECT_GT(e, 0.0);
    EXPECT_LT(e, 2.0);
}

TEST_F(TrefftzGradientTest, SpanEfficiencyZeroDrag) {
    double e = compute_span_efficiency(0.5, 0.0, 0.0, ref_area_, ref_span_);
    EXPECT_DOUBLE_EQ(e, 0.0);
}

TEST_F(TrefftzGradientTest, TrefftzVelocityPoint) {
    double vy = 0.0, vz = 0.0;
    compute_trefftz_velocity_point(1.0, 0.0, 0.0, 0.5, 0.0, 0.0, &vy, &vz);
    EXPECT_TRUE(std::isfinite(vy));
    EXPECT_TRUE(std::isfinite(vz));
    EXPECT_NE(vz, 0.0);
}

TEST_F(TrefftzGradientTest, CLGradientMatchesFD) {
    std::vector<double> d_gamma(num_strips_, 0.0);
    analytical_trefftz_circulation_gradient(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, 0, d_gamma.data());

    for (int k = 0; k < num_strips_; ++k) {
        auto gamma_p = gamma_;
        auto gamma_m = gamma_;
        gamma_p[k] += kFDStep;
        gamma_m[k] -= kFDStep;

        double CL_p, CY_p, CDi_p, CL_m, CY_m, CDi_m;
        compute_trefftz_drag_kernel(
            gamma_p.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_p, &CY_p, &CDi_p);
        compute_trefftz_drag_kernel(
            gamma_m.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_m, &CY_m, &CDi_m);

        double fd = (CL_p - CL_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_gamma[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Strip " << k;
    }
}

TEST_F(TrefftzGradientTest, CYGradientMatchesFD) {
    std::vector<double> d_gamma(num_strips_, 0.0);
    analytical_trefftz_circulation_gradient(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, 1, d_gamma.data());

    for (int k = 0; k < num_strips_; ++k) {
        auto gamma_p = gamma_;
        auto gamma_m = gamma_;
        gamma_p[k] += kFDStep;
        gamma_m[k] -= kFDStep;

        double CL_p, CY_p, CDi_p, CL_m, CY_m, CDi_m;
        compute_trefftz_drag_kernel(
            gamma_p.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_p, &CY_p, &CDi_p);
        compute_trefftz_drag_kernel(
            gamma_m.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_m, &CY_m, &CDi_m);

        double fd = (CY_p - CY_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_gamma[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Strip " << k;
    }
}

TEST_F(TrefftzGradientTest, CDiGradientMatchesFD) {
    std::vector<double> d_gamma(num_strips_, 0.0);
    analytical_trefftz_circulation_gradient(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, 2, d_gamma.data());

    for (int k = 0; k < num_strips_; ++k) {
        auto gamma_p = gamma_;
        auto gamma_m = gamma_;
        gamma_p[k] += kFDStep;
        gamma_m[k] -= kFDStep;

        double CL_p, CY_p, CDi_p, CL_m, CY_m, CDi_m;
        compute_trefftz_drag_kernel(
            gamma_p.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_p, &CY_p, &CDi_p);
        compute_trefftz_drag_kernel(
            gamma_m.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_m, &CY_m, &CDi_m);

        double fd = (CDi_p - CDi_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_gamma[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Strip " << k;
    }
}

TEST_F(TrefftzGradientTest, MultiSurfaceCoreRadius) {
    // Two surfaces with different indices -> non-zero core radius
    surface_idx_ = {0, 0, 1};

    std::vector<double> d_gamma(num_strips_, 0.0);
    analytical_trefftz_circulation_gradient(
        gamma_.data(), te_left_.data(), te_right_.data(), chord_.data(),
        surface_idx_.data(), num_strips_, ref_area_, ref_span_,
        core_chord_frac_, core_width_frac_, 2, d_gamma.data());

    // FD validation still holds
    for (int k = 0; k < num_strips_; ++k) {
        auto gamma_p = gamma_;
        auto gamma_m = gamma_;
        gamma_p[k] += kFDStep;
        gamma_m[k] -= kFDStep;

        double CL_p, CY_p, CDi_p, CL_m, CY_m, CDi_m;
        compute_trefftz_drag_kernel(
            gamma_p.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_p, &CY_p, &CDi_p);
        compute_trefftz_drag_kernel(
            gamma_m.data(), te_left_.data(), te_right_.data(), chord_.data(),
            surface_idx_.data(), num_strips_, ref_area_, ref_span_,
            core_chord_frac_, core_width_frac_, &CL_m, &CY_m, &CDi_m);

        double fd = (CDi_p - CDi_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_gamma[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Strip " << k;
    }
}

}  // namespace
}  // namespace pegasus::gradients
