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

#include "pegasus/gradients/analytical_biot_savart.h"

#include <cmath>

#include <gtest/gtest.h>

namespace pegasus::gradients {
namespace {

constexpr double kFDStep = 1e-6;
constexpr double kFDTol = 1e-4;

// Helper: perturb one component, compute velocity, return finite-difference VJP
void fd_vjp_field_point(const double* fp, const double* bl, const double* br,
                        double pg, double cr, bool ib, const double* seed,
                        double* result) {
    for (int k = 0; k < 3; ++k) {
        double fp_p[3] = {fp[0], fp[1], fp[2]};
        double fp_m[3] = {fp[0], fp[1], fp[2]};
        fp_p[k] += kFDStep;
        fp_m[k] -= kFDStep;

        double v_p[3], v_m[3];
        compute_horseshoe_velocity_kernel(fp_p, bl, br, pg, cr, ib, v_p);
        compute_horseshoe_velocity_kernel(fp_m, bl, br, pg, cr, ib, v_m);

        result[k] = 0.0;
        for (int m = 0; m < 3; ++m) {
            result[k] += seed[m] * (v_p[m] - v_m[m]) / (2.0 * kFDStep);
        }
    }
}

void fd_vjp_bound_left(const double* fp, const double* bl, const double* br,
                       double pg, double cr, bool ib, const double* seed,
                       double* result) {
    for (int k = 0; k < 3; ++k) {
        double bl_p[3] = {bl[0], bl[1], bl[2]};
        double bl_m[3] = {bl[0], bl[1], bl[2]};
        bl_p[k] += kFDStep;
        bl_m[k] -= kFDStep;

        double v_p[3], v_m[3];
        compute_horseshoe_velocity_kernel(fp, bl_p, br, pg, cr, ib, v_p);
        compute_horseshoe_velocity_kernel(fp, bl_m, br, pg, cr, ib, v_m);

        result[k] = 0.0;
        for (int m = 0; m < 3; ++m) {
            result[k] += seed[m] * (v_p[m] - v_m[m]) / (2.0 * kFDStep);
        }
    }
}

void fd_vjp_bound_right(const double* fp, const double* bl, const double* br,
                        double pg, double cr, bool ib, const double* seed,
                        double* result) {
    for (int k = 0; k < 3; ++k) {
        double br_p[3] = {br[0], br[1], br[2]};
        double br_m[3] = {br[0], br[1], br[2]};
        br_p[k] += kFDStep;
        br_m[k] -= kFDStep;

        double v_p[3], v_m[3];
        compute_horseshoe_velocity_kernel(fp, bl, br_p, pg, cr, ib, v_p);
        compute_horseshoe_velocity_kernel(fp, bl, br_m, pg, cr, ib, v_m);

        result[k] = 0.0;
        for (int m = 0; m < 3; ++m) {
            result[k] += seed[m] * (v_p[m] - v_m[m]) / (2.0 * kFDStep);
        }
    }
}

// Typical horseshoe geometry: field point above bound vortex
class BiotSavartVJPTest : public ::testing::Test {
protected:
    double fp_[3] = {0.5, 0.5, 0.0};
    double bl_[3] = {0.25, 0.0, 0.0};
    double br_[3] = {0.25, 1.0, 0.0};
    double seed_[3] = {0.0, 0.0, 1.0};
};

TEST_F(BiotSavartVJPTest, PrimalNonZero) {
    double vel[3];
    compute_horseshoe_velocity_kernel(fp_, bl_, br_, 1.0, 0.01, true, vel);
    EXPECT_TRUE(std::isfinite(vel[0]));
    EXPECT_TRUE(std::isfinite(vel[1]));
    EXPECT_TRUE(std::isfinite(vel[2]));
    // At least one component should be non-zero
    double mag = std::sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    EXPECT_GT(mag, 1e-10);
}

TEST_F(BiotSavartVJPTest, FieldPointVJPMatchesFD) {
    double analytical[3], fd[3];
    horseshoe_velocity_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, true, seed_,
                                       analytical);
    fd_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, true, seed_, fd);

    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8))
            << "Component " << k;
    }
}

TEST_F(BiotSavartVJPTest, BoundLeftVJPMatchesFD) {
    double analytical[3], fd[3];
    horseshoe_velocity_vjp_bound_left(fp_, bl_, br_, 1.0, 0.01, true, seed_,
                                      analytical);
    fd_vjp_bound_left(fp_, bl_, br_, 1.0, 0.01, true, seed_, fd);

    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8))
            << "Component " << k;
    }
}

TEST_F(BiotSavartVJPTest, BoundRightVJPMatchesFD) {
    double analytical[3], fd[3];
    horseshoe_velocity_vjp_bound_right(fp_, bl_, br_, 1.0, 0.01, true, seed_,
                                       analytical);
    fd_vjp_bound_right(fp_, bl_, br_, 1.0, 0.01, true, seed_, fd);

    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8))
            << "Component " << k;
    }
}

TEST_F(BiotSavartVJPTest, WithoutBoundVortex) {
    double analytical[3], fd[3];

    horseshoe_velocity_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, false, seed_,
                                       analytical);
    fd_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, false, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8));
    }

    horseshoe_velocity_vjp_bound_left(fp_, bl_, br_, 1.0, 0.01, false, seed_,
                                      analytical);
    fd_vjp_bound_left(fp_, bl_, br_, 1.0, 0.01, false, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8));
    }

    horseshoe_velocity_vjp_bound_right(fp_, bl_, br_, 1.0, 0.01, false, seed_,
                                       analytical);
    fd_vjp_bound_right(fp_, bl_, br_, 1.0, 0.01, false, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8));
    }
}

TEST_F(BiotSavartVJPTest, PrandtlGlauertScaling) {
    double pg = 0.8;  // Mach ~0.6
    double analytical[3], fd[3];

    horseshoe_velocity_vjp_field_point(fp_, bl_, br_, pg, 0.01, true, seed_,
                                       analytical);
    fd_vjp_field_point(fp_, bl_, br_, pg, 0.01, true, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-4));
    }

    horseshoe_velocity_vjp_bound_left(fp_, bl_, br_, pg, 0.01, true, seed_,
                                      analytical);
    fd_vjp_bound_left(fp_, bl_, br_, pg, 0.01, true, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-4));
    }

    horseshoe_velocity_vjp_bound_right(fp_, bl_, br_, pg, 0.01, true, seed_,
                                       analytical);
    fd_vjp_bound_right(fp_, bl_, br_, pg, 0.01, true, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-4));
    }
}

TEST_F(BiotSavartVJPTest, LargerCoreRadius) {
    double cr = 0.1;
    double analytical[3], fd[3];

    horseshoe_velocity_vjp_field_point(fp_, bl_, br_, 1.0, cr, true, seed_,
                                       analytical);
    fd_vjp_field_point(fp_, bl_, br_, 1.0, cr, true, seed_, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8));
    }
}

TEST(BiotSavartAICTest, AICElementKernel) {
    double fp[3] = {0.75, 0.5, 0.0};
    double bl[3] = {0.25, 0.0, 0.0};
    double br[3] = {0.25, 1.0, 0.0};
    double normal[3] = {0.0, 0.0, 1.0};

    double aic =
        compute_aic_element_kernel(fp, normal, bl, br, 1.0, 0.01, true);

    // AIC should be finite and non-zero
    EXPECT_TRUE(std::isfinite(aic));
    EXPECT_NE(aic, 0.0);

    // Verify it equals V . normal
    double vel[3];
    compute_horseshoe_velocity_kernel(fp, bl, br, 1.0, 0.01, true, vel);
    double dot = vel[0] * normal[0] + vel[1] * normal[1] + vel[2] * normal[2];
    EXPECT_DOUBLE_EQ(aic, dot);
}

TEST(BiotSavartAICTest, AICElementKernelWithoutBound) {
    double fp[3] = {0.75, 0.5, 0.0};
    double bl[3] = {0.25, 0.0, 0.0};
    double br[3] = {0.25, 1.0, 0.0};
    double normal[3] = {0.0, 0.0, 1.0};

    double aic_with =
        compute_aic_element_kernel(fp, normal, bl, br, 1.0, 0.01, true);
    double aic_without =
        compute_aic_element_kernel(fp, normal, bl, br, 1.0, 0.01, false);

    // Without bound vortex should differ
    EXPECT_NE(aic_with, aic_without);
    EXPECT_TRUE(std::isfinite(aic_without));
}

TEST_F(BiotSavartVJPTest, ArbitrarySeedDirection) {
    double seed[3] = {1.0, -0.5, 0.3};
    double analytical[3], fd[3];

    horseshoe_velocity_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, true, seed,
                                       analytical);
    fd_vjp_field_point(fp_, bl_, br_, 1.0, 0.01, true, seed, fd);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(analytical[k], fd[k], kFDTol * (std::abs(fd[k]) + 1e-8));
    }
}

}  // namespace
}  // namespace pegasus::gradients
