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

#include "pegasus/gradients/analytical_aic.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace pegasus::gradients {
namespace {

constexpr double kFDStep = 1e-6;
constexpr double kFDTol = 1e-4;

// ============================================================================
// AIC Element Gradient Tests
// ============================================================================

class AICElementGradientTest : public ::testing::Test {
protected:
    double rv1_[3] = {0.25, 0.0, 0.0};
    double rv2_[3] = {0.25, 1.0, 0.0};
    double rc_[3] = {0.75, 0.5, 0.0};
    double normal_[3] = {0.0, 0.0, 1.0};
    double core_radius_ = 0.01;
    double pg_ = 1.0;
    bool include_bound_ = true;
};

TEST_F(AICElementGradientTest, GradientsMatchFD) {
    double d_rv1[3], d_rv2[3], d_rc[3], d_normal[3];
    analytical_aic_element_gradient(rv1_, rv2_, rc_, normal_, core_radius_, pg_,
                                    include_bound_, d_rv1, d_rv2, d_rc,
                                    d_normal);

    // FD for rc
    for (int k = 0; k < 3; ++k) {
        double rc_p[3] = {rc_[0], rc_[1], rc_[2]};
        double rc_m[3] = {rc_[0], rc_[1], rc_[2]};
        rc_p[k] += kFDStep;
        rc_m[k] -= kFDStep;
        double aic_p = compute_aic_element_kernel(
            rc_p, normal_, rv1_, rv2_, pg_, core_radius_, include_bound_);
        double aic_m = compute_aic_element_kernel(
            rc_m, normal_, rv1_, rv2_, pg_, core_radius_, include_bound_);
        double fd = (aic_p - aic_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rc[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "d_rc[" << k << "]";
    }

    // FD for rv1
    for (int k = 0; k < 3; ++k) {
        double rv1_p[3] = {rv1_[0], rv1_[1], rv1_[2]};
        double rv1_m[3] = {rv1_[0], rv1_[1], rv1_[2]};
        rv1_p[k] += kFDStep;
        rv1_m[k] -= kFDStep;
        double aic_p = compute_aic_element_kernel(
            rc_, normal_, rv1_p, rv2_, pg_, core_radius_, include_bound_);
        double aic_m = compute_aic_element_kernel(
            rc_, normal_, rv1_m, rv2_, pg_, core_radius_, include_bound_);
        double fd = (aic_p - aic_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rv1[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "d_rv1[" << k << "]";
    }

    // FD for rv2
    for (int k = 0; k < 3; ++k) {
        double rv2_p[3] = {rv2_[0], rv2_[1], rv2_[2]};
        double rv2_m[3] = {rv2_[0], rv2_[1], rv2_[2]};
        rv2_p[k] += kFDStep;
        rv2_m[k] -= kFDStep;
        double aic_p = compute_aic_element_kernel(
            rc_, normal_, rv1_, rv2_p, pg_, core_radius_, include_bound_);
        double aic_m = compute_aic_element_kernel(
            rc_, normal_, rv1_, rv2_m, pg_, core_radius_, include_bound_);
        double fd = (aic_p - aic_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rv2[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "d_rv2[" << k << "]";
    }

    // d_normal should equal velocity
    double vel[3];
    compute_horseshoe_velocity_kernel(rc_, rv1_, rv2_, pg_, core_radius_,
                                      include_bound_, vel);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(d_normal[k], vel[k], 1e-12) << "d_normal[" << k << "]";
    }
}

TEST_F(AICElementGradientTest, WithoutBound) {
    include_bound_ = false;
    double d_rv1[3], d_rv2[3], d_rc[3], d_normal[3];
    analytical_aic_element_gradient(rv1_, rv2_, rc_, normal_, core_radius_, pg_,
                                    include_bound_, d_rv1, d_rv2, d_rc,
                                    d_normal);

    // FD for rc
    for (int k = 0; k < 3; ++k) {
        double rc_p[3] = {rc_[0], rc_[1], rc_[2]};
        double rc_m[3] = {rc_[0], rc_[1], rc_[2]};
        rc_p[k] += kFDStep;
        rc_m[k] -= kFDStep;
        double aic_p = compute_aic_element_kernel(rc_p, normal_, rv1_, rv2_,
                                                  pg_, core_radius_, false);
        double aic_m = compute_aic_element_kernel(rc_m, normal_, rv1_, rv2_,
                                                  pg_, core_radius_, false);
        double fd = (aic_p - aic_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rc[k], fd, kFDTol * (std::abs(fd) + 1e-8));
    }
}

TEST_F(AICElementGradientTest, PrandtlGlauert) {
    pg_ = 0.8;
    double d_rv1[3], d_rv2[3], d_rc[3], d_normal[3];
    analytical_aic_element_gradient(rv1_, rv2_, rc_, normal_, core_radius_, pg_,
                                    include_bound_, d_rv1, d_rv2, d_rc,
                                    d_normal);

    for (int k = 0; k < 3; ++k) {
        double rc_p[3] = {rc_[0], rc_[1], rc_[2]};
        double rc_m[3] = {rc_[0], rc_[1], rc_[2]};
        rc_p[k] += kFDStep;
        rc_m[k] -= kFDStep;
        double aic_p = compute_aic_element_kernel(rc_p, normal_, rv1_, rv2_,
                                                  pg_, core_radius_, true);
        double aic_m = compute_aic_element_kernel(rc_m, normal_, rv1_, rv2_,
                                                  pg_, core_radius_, true);
        double fd = (aic_p - aic_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rc[k], fd, kFDTol * (std::abs(fd) + 1e-8));
    }
}

// ============================================================================
// AIC Row/Column Gradient Tests
// ============================================================================

class AICRowColumnTest : public ::testing::Test {
protected:
    void SetUp() override {
        num_vortices_ = 2;
        // Two horseshoe vortices side by side
        rv1_ = {0.25, 0.0, 0.0, 0.25, 0.5, 0.0};
        rv2_ = {0.25, 0.5, 0.0, 0.25, 1.0, 0.0};
        rc_ = {0.75, 0.25, 0.0, 0.75, 0.75, 0.0};
        normal_ = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
        chord_ = {1.0, 1.0};
        component_idx_ = {0, 0};
        weights_ = {1.0, 0.5};
    }

    int num_vortices_;
    std::vector<double> rv1_, rv2_, rc_, normal_, chord_;
    std::vector<int> component_idx_;
    std::vector<double> weights_;
    double core_chord_frac_ = 0.25;
    double core_width_frac_ = 0.05;
    double beta_ = 1.0;
};

TEST_F(AICRowColumnTest, RowGradientMatchesFD) {
    int row_idx = 0;
    int control_comp = component_idx_[row_idx];

    double d_rc[3], d_normal[3];
    analytical_aic_row_gradient_rc(
        rv1_.data(), rv2_.data(), &rc_[3 * row_idx], &normal_[3 * row_idx],
        chord_.data(), component_idx_.data(), num_vortices_, row_idx,
        control_comp, core_chord_frac_, core_width_frac_, beta_,
        weights_.data(), d_rc, d_normal);

    // FD for rc
    for (int k = 0; k < 3; ++k) {
        auto rc_pert = rc_;
        rc_pert[3 * row_idx + k] += kFDStep;

        // Compute sum_j AIC[i,j]*weights[j] with perturbed rc
        double val_p = 0.0;
        for (int j = 0; j < num_vortices_; ++j) {
            bool ib = (j != row_idx);
            double sw =
                std::sqrt(std::pow(rv2_[3 * j + 1] - rv1_[3 * j + 1], 2) +
                          std::pow(rv2_[3 * j + 2] - rv1_[3 * j + 2], 2));
            double cr = 0.0001 * sw;  // same component
            val_p +=
                weights_[j] * compute_aic_element_kernel(
                                  &rc_pert[3 * row_idx], &normal_[3 * row_idx],
                                  &rv1_[3 * j], &rv2_[3 * j], beta_, cr, ib);
        }

        rc_pert[3 * row_idx + k] -= 2.0 * kFDStep;
        double val_m = 0.0;
        for (int j = 0; j < num_vortices_; ++j) {
            bool ib = (j != row_idx);
            double sw =
                std::sqrt(std::pow(rv2_[3 * j + 1] - rv1_[3 * j + 1], 2) +
                          std::pow(rv2_[3 * j + 2] - rv1_[3 * j + 2], 2));
            double cr = 0.0001 * sw;
            val_m +=
                weights_[j] * compute_aic_element_kernel(
                                  &rc_pert[3 * row_idx], &normal_[3 * row_idx],
                                  &rv1_[3 * j], &rv2_[3 * j], beta_, cr, ib);
        }

        double fd = (val_p - val_m) / (2.0 * kFDStep);
        EXPECT_NEAR(d_rc[k], fd, kFDTol * (std::abs(fd) + 1e-8))
            << "d_rc[" << k << "]";
    }
}

TEST_F(AICRowColumnTest, ColumnGradientFinite) {
    double d_rv1[3], d_rv2[3];
    analytical_aic_column_gradient_rv(
        rv1_.data(), rv2_.data(), rc_.data(), normal_.data(), chord_.data(),
        component_idx_.data(), num_vortices_, 0, core_chord_frac_,
        core_width_frac_, beta_, weights_.data(), d_rv1, d_rv2);

    for (int k = 0; k < 3; ++k) {
        EXPECT_TRUE(std::isfinite(d_rv1[k]));
        EXPECT_TRUE(std::isfinite(d_rv2[k]));
    }
}

TEST_F(AICRowColumnTest, DifferentComponents) {
    // Put vortices on different components
    component_idx_ = {0, 1};

    double d_rc[3], d_normal[3];
    analytical_aic_row_gradient_rc(
        rv1_.data(), rv2_.data(), &rc_[0], &normal_[0], chord_.data(),
        component_idx_.data(), num_vortices_, 0, 0, core_chord_frac_,
        core_width_frac_, beta_, weights_.data(), d_rc, d_normal);

    for (int k = 0; k < 3; ++k) {
        EXPECT_TRUE(std::isfinite(d_rc[k]));
        EXPECT_TRUE(std::isfinite(d_normal[k]));
    }
}

TEST_F(AICRowColumnTest, VelocityRowGradientFinite) {
    double dJ_dV[3] = {0.1, -0.2, 0.5};
    double d_rc[3];
    analytical_velocity_row_gradient_rc(
        rv1_.data(), rv2_.data(), &rc_[0], chord_.data(), component_idx_.data(),
        num_vortices_, 0, 0, core_chord_frac_, core_width_frac_, beta_,
        weights_.data(), dJ_dV, d_rc);

    for (int k = 0; k < 3; ++k) {
        EXPECT_TRUE(std::isfinite(d_rc[k]));
    }
}

TEST_F(AICRowColumnTest, VelocityColumnGradientFinite) {
    std::vector<double> dJ_dV(3 * num_vortices_, 0.0);
    dJ_dV[2] = 1.0;
    dJ_dV[5] = 0.5;

    double d_rv1[3], d_rv2[3];
    analytical_velocity_column_gradient_rv(
        rv1_.data(), rv2_.data(), rc_.data(), chord_.data(),
        component_idx_.data(), num_vortices_, 0, core_chord_frac_,
        core_width_frac_, beta_, 1.0, dJ_dV.data(), d_rv1, d_rv2);

    for (int k = 0; k < 3; ++k) {
        EXPECT_TRUE(std::isfinite(d_rv1[k]));
        EXPECT_TRUE(std::isfinite(d_rv2[k]));
    }
}

TEST_F(AICRowColumnTest, VelocityColumnGradientZeroGamma) {
    std::vector<double> dJ_dV(3 * num_vortices_, 1.0);
    double d_rv1[3], d_rv2[3];

    analytical_velocity_column_gradient_rv(
        rv1_.data(), rv2_.data(), rc_.data(), chord_.data(),
        component_idx_.data(), num_vortices_, 0, core_chord_frac_,
        core_width_frac_, beta_, 0.0, dJ_dV.data(), d_rv1, d_rv2);

    for (int k = 0; k < 3; ++k) {
        EXPECT_DOUBLE_EQ(d_rv1[k], 0.0);
        EXPECT_DOUBLE_EQ(d_rv2[k], 0.0);
    }
}

}  // namespace
}  // namespace pegasus::gradients
