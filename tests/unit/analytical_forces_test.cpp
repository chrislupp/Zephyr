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

#include "pegasus/gradients/analytical_forces.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace pegasus::gradients {
namespace {

constexpr double kFDStep = 1e-7;
constexpr double kFDTol = 1e-4;

class ForceGradientTest : public ::testing::Test {
protected:
    void SetUp() override {
        num_vortices_ = 2;
        alpha_ = 5.0 * M_PI / 180.0;
        beta_ = 0.0;
        ref_area_ = 10.0;

        circulation_ = {0.1, 0.08};
        induced_velocity_ = {0.0, 0.0, -0.005, 0.0, 0.0, -0.004};
        bound_left_ = {0.25, 0.0, 0.0, 0.25, 0.5, 0.0};
        bound_right_ = {0.25, 0.5, 0.0, 0.25, 1.0, 0.0};
        control_points_ = {0.75, 0.25, 0.0, 0.75, 0.75, 0.0};
        moment_arm_ = {0.5, 0.25, 0.0, 0.5, 0.75, 0.0};
    }

    int num_vortices_;
    double alpha_, beta_, ref_area_;
    std::vector<double> circulation_, induced_velocity_;
    std::vector<double> bound_left_, bound_right_;
    std::vector<double> control_points_, moment_arm_;
};

TEST_F(ForceGradientTest, PrimalKernelComputesAllOutputs) {
    double result[6];
    compute_force_coefficients_kernel(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result);

    // CL should be positive at positive alpha
    EXPECT_GT(result[0], 0.0);  // CL > 0
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(std::isfinite(result[i])) << "Output " << i;
    }
}

TEST_F(ForceGradientTest, FreestreamDirection) {
    double v[3];
    compute_freestream_direction(0.0, 0.0, v);
    EXPECT_NEAR(v[0], 1.0, 1e-12);
    EXPECT_NEAR(v[1], 0.0, 1e-12);
    EXPECT_NEAR(v[2], 0.0, 1e-12);

    compute_freestream_direction(M_PI / 4.0, 0.0, v);
    EXPECT_NEAR(v[0], std::cos(M_PI / 4.0), 1e-12);
    EXPECT_NEAR(v[2], std::sin(M_PI / 4.0), 1e-12);
}

TEST_F(ForceGradientTest, FreestreamAlphaGradient) {
    double dv[3];
    compute_freestream_alpha_gradient(alpha_, beta_, dv);

    // FD check
    double v_p[3], v_m[3];
    compute_freestream_direction(alpha_ + kFDStep, beta_, v_p);
    compute_freestream_direction(alpha_ - kFDStep, beta_, v_m);
    for (int k = 0; k < 3; ++k) {
        double fd = (v_p[k] - v_m[k]) / (2.0 * kFDStep);
        EXPECT_NEAR(dv[k], fd, 1e-8);
    }
}

TEST_F(ForceGradientTest, FreestreamBetaGradient) {
    beta_ = 3.0 * M_PI / 180.0;
    double dv[3];
    compute_freestream_beta_gradient(alpha_, beta_, dv);

    double v_p[3], v_m[3];
    compute_freestream_direction(alpha_, beta_ + kFDStep, v_p);
    compute_freestream_direction(alpha_, beta_ - kFDStep, v_m);
    for (int k = 0; k < 3; ++k) {
        double fd = (v_p[k] - v_m[k]) / (2.0 * kFDStep);
        EXPECT_NEAR(dv[k], fd, 1e-8);
    }
}

TEST_F(ForceGradientTest, BodyToStability) {
    double CL, CD;
    body_to_stability(1.0, 0.0, 0.0, &CL, &CD);
    EXPECT_NEAR(CL, 0.0, 1e-12);
    EXPECT_NEAR(CD, 1.0, 1e-12);

    body_to_stability(0.0, 1.0, 0.0, &CL, &CD);
    EXPECT_NEAR(CL, 1.0, 1e-12);
    EXPECT_NEAR(CD, 0.0, 1e-12);
}

TEST_F(ForceGradientTest, AlphaGradientMatchesFD) {
    double d_result[6];
    analytical_force_alpha_gradient(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, d_result);

    double result_p[6], result_m[6];
    compute_force_coefficients_kernel(
        alpha_ + kFDStep, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_p);
    compute_force_coefficients_kernel(
        alpha_ - kFDStep, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_m);

    for (int i = 0; i < 6; ++i) {
        double fd = (result_p[i] - result_m[i]) / (2.0 * kFDStep);
        EXPECT_NEAR(d_result[i], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Output " << i;
    }
}

TEST_F(ForceGradientTest, BetaGradientMatchesFD) {
    beta_ = 3.0 * M_PI / 180.0;
    double d_result[6];
    analytical_force_beta_gradient(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, d_result);

    double result_p[6], result_m[6];
    compute_force_coefficients_kernel(
        alpha_, beta_ + kFDStep, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_p);
    compute_force_coefficients_kernel(
        alpha_, beta_ - kFDStep, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_m);

    for (int i = 0; i < 6; ++i) {
        double fd = (result_p[i] - result_m[i]) / (2.0 * kFDStep);
        EXPECT_NEAR(d_result[i], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Output " << i;
    }
}

TEST_F(ForceGradientTest, CirculationGradientMatchesFD) {
    for (int output_idx = 0; output_idx < 6; ++output_idx) {
        std::vector<double> d_circ(num_vortices_, 0.0);
        analytical_force_circulation_gradient(
            alpha_, beta_, circulation_.data(), induced_velocity_.data(),
            bound_left_.data(), bound_right_.data(), control_points_.data(),
            moment_arm_.data(), num_vortices_, ref_area_, output_idx,
            d_circ.data());

        for (int j = 0; j < num_vortices_; ++j) {
            auto circ_p = circulation_;
            auto circ_m = circulation_;
            circ_p[j] += kFDStep;
            circ_m[j] -= kFDStep;

            double result_p[6], result_m[6];
            compute_force_coefficients_kernel(
                alpha_, beta_, circ_p.data(), induced_velocity_.data(),
                bound_left_.data(), bound_right_.data(), control_points_.data(),
                moment_arm_.data(), num_vortices_, ref_area_, result_p);
            compute_force_coefficients_kernel(
                alpha_, beta_, circ_m.data(), induced_velocity_.data(),
                bound_left_.data(), bound_right_.data(), control_points_.data(),
                moment_arm_.data(), num_vortices_, ref_area_, result_m);

            double fd =
                (result_p[output_idx] - result_m[output_idx]) / (2.0 * kFDStep);
            EXPECT_NEAR(d_circ[j], fd, kFDTol * (std::abs(fd) + 1e-8))
                << "output=" << output_idx << " vortex=" << j;
        }
    }
}

TEST_F(ForceGradientTest, InducedVelocityGradientMatchesFD) {
    for (int output_idx = 0; output_idx < 6; ++output_idx) {
        std::vector<double> d_vi(3 * num_vortices_, 0.0);
        analytical_force_induced_velocity_gradient(
            alpha_, beta_, circulation_.data(), induced_velocity_.data(),
            bound_left_.data(), bound_right_.data(), control_points_.data(),
            moment_arm_.data(), num_vortices_, ref_area_, output_idx,
            d_vi.data());

        for (int j = 0; j < 3 * num_vortices_; ++j) {
            auto vi_p = induced_velocity_;
            auto vi_m = induced_velocity_;
            vi_p[j] += kFDStep;
            vi_m[j] -= kFDStep;

            double result_p[6], result_m[6];
            compute_force_coefficients_kernel(
                alpha_, beta_, circulation_.data(), vi_p.data(),
                bound_left_.data(), bound_right_.data(), control_points_.data(),
                moment_arm_.data(), num_vortices_, ref_area_, result_p);
            compute_force_coefficients_kernel(
                alpha_, beta_, circulation_.data(), vi_m.data(),
                bound_left_.data(), bound_right_.data(), control_points_.data(),
                moment_arm_.data(), num_vortices_, ref_area_, result_m);

            double fd =
                (result_p[output_idx] - result_m[output_idx]) / (2.0 * kFDStep);
            EXPECT_NEAR(d_vi[j], fd, kFDTol * (std::abs(fd) + 1e-8))
                << "output=" << output_idx << " vi_idx=" << j;
        }
    }
}

TEST_F(ForceGradientTest, GeometryGradientMatchesFD) {
    // Forward-mode geometry gradient: perturb bound_left, bound_right,
    // moment_arm
    std::vector<double> d_bl(3 * num_vortices_, 0.0);
    std::vector<double> d_br(3 * num_vortices_, 0.0);
    std::vector<double> d_ma(3 * num_vortices_, 0.0);

    // Perturb first bound left y component
    d_bl[1] = 1.0;

    double d_result[6];
    analytical_force_geometry_gradient(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), moment_arm_.data(),
        num_vortices_, ref_area_, d_bl.data(), d_br.data(), d_ma.data(),
        d_result);

    // FD check
    auto bl_p = bound_left_;
    auto bl_m = bound_left_;
    bl_p[1] += kFDStep;
    bl_m[1] -= kFDStep;

    double result_p[6], result_m[6];
    compute_force_coefficients_kernel(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bl_p.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_p);
    compute_force_coefficients_kernel(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bl_m.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, result_m);

    for (int i = 0; i < 6; ++i) {
        double fd = (result_p[i] - result_m[i]) / (2.0 * kFDStep);
        EXPECT_NEAR(d_result[i], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "Output " << i;
    }
}

TEST_F(ForceGradientTest, ZeroAlpha) {
    alpha_ = 0.0;
    double d_result[6];
    analytical_force_alpha_gradient(
        alpha_, beta_, circulation_.data(), induced_velocity_.data(),
        bound_left_.data(), bound_right_.data(), control_points_.data(),
        moment_arm_.data(), num_vortices_, ref_area_, d_result);

    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(std::isfinite(d_result[i])) << "Output " << i;
    }
}

}  // namespace
}  // namespace pegasus::gradients
