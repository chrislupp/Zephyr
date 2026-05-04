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

#include "pegasus/aero/aerodynamic_results.h"

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTol = 1.0e-10;

// ============================================================================
// AerodynamicResults Tests
// ============================================================================

TEST(AerodynamicResultsTest, DefaultState) {
    AerodynamicResults r;
    EXPECT_FALSE(r.is_valid);
    EXPECT_NEAR(r.CL, 0.0, kTol);
    EXPECT_NEAR(r.CD, 0.0, kTol);
    EXPECT_NEAR(r.reference_area, 1.0, kTol);
    EXPECT_NEAR(r.reference_chord, 1.0, kTol);
    EXPECT_NEAR(r.reference_span, 1.0, kTol);
}

TEST(AerodynamicResultsTest, GetAspectRatio) {
    AerodynamicResults r;
    r.reference_span = 10.0;
    r.reference_area = 20.0;
    // AR = b^2 / S = 100 / 20 = 5
    EXPECT_NEAR(r.GetAspectRatio(), 5.0, kTol);
}

TEST(AerodynamicResultsTest, GetAspectRatioSquareWing) {
    AerodynamicResults r;
    r.reference_span = 4.0;
    r.reference_area = 4.0;
    EXPECT_NEAR(r.GetAspectRatio(), 4.0, kTol);
}

TEST(AerodynamicResultsTest, GetIdealInducedDrag) {
    AerodynamicResults r;
    r.reference_span = 10.0;
    r.reference_area = 10.0;
    r.CL = 0.5;
    // AR = 100/10 = 10, CDi_ideal = CL^2 / (pi * AR) = 0.25 / (pi*10)
    Scalar expected = 0.25 / (kPi * 10.0);
    EXPECT_NEAR(r.GetIdealInducedDrag(), expected, 1e-8);
}

TEST(AerodynamicResultsTest, GetIdealInducedDragZeroCL) {
    AerodynamicResults r;
    r.CL = 0.0;
    EXPECT_NEAR(r.GetIdealInducedDrag(), 0.0, kTol);
}

TEST(AerodynamicResultsTest, Clear) {
    AerodynamicResults r;
    r.is_valid = true;
    r.CL = 0.5;
    r.CD = 0.01;
    r.CX = 0.1;
    r.CY = 0.02;
    r.CZ = -0.5;
    r.Cl = 0.001;
    r.Cm = -0.05;
    r.Cn = 0.002;
    r.CDi_nearfield = 0.008;
    r.CDv = 0.003;
    r.CL_ff = 0.49;
    r.CY_ff = 0.01;
    r.CDi_ff = 0.007;
    r.span_efficiency = 0.95;
    r.CF = Vec3(1, 2, 3);
    r.CM = Vec3(4, 5, 6);

    SurfaceLoads loads;
    loads.name = "wing";
    r.surface_loads.push_back(loads);

    r.Clear();

    EXPECT_FALSE(r.is_valid);
    EXPECT_NEAR(r.CL, 0.0, kTol);
    EXPECT_NEAR(r.CD, 0.0, kTol);
    EXPECT_NEAR(r.CX, 0.0, kTol);
    EXPECT_NEAR(r.CY, 0.0, kTol);
    EXPECT_NEAR(r.CZ, 0.0, kTol);
    EXPECT_NEAR(r.Cl, 0.0, kTol);
    EXPECT_NEAR(r.Cm, 0.0, kTol);
    EXPECT_NEAR(r.Cn, 0.0, kTol);
    EXPECT_NEAR(r.CDi_nearfield, 0.0, kTol);
    EXPECT_NEAR(r.CDv, 0.0, kTol);
    EXPECT_NEAR(r.CL_ff, 0.0, kTol);
    EXPECT_NEAR(r.CY_ff, 0.0, kTol);
    EXPECT_NEAR(r.CDi_ff, 0.0, kTol);
    EXPECT_NEAR(r.span_efficiency, 0.0, kTol);
    EXPECT_NEAR(r.CF.norm(), 0.0, kTol);
    EXPECT_NEAR(r.CM.norm(), 0.0, kTol);
    EXPECT_TRUE(r.surface_loads.empty());
}

// ============================================================================
// StripLoads Tests
// ============================================================================

TEST(StripLoadsTest, DefaultState) {
    StripLoads s;
    EXPECT_EQ(s.strip_index, 0);
    EXPECT_NEAR(s.chord, 0.0, kTol);
    EXPECT_NEAR(s.width, 0.0, kTol);
    EXPECT_NEAR(s.area, 0.0, kTol);
    EXPECT_NEAR(s.cl, 0.0, kTol);
    EXPECT_NEAR(s.cd, 0.0, kTol);
    EXPECT_NEAR(s.cy, 0.0, kTol);
    EXPECT_NEAR(s.cnc, 0.0, kTol);
    EXPECT_NEAR(s.alpha_induced, 0.0, kTol);
}

// ============================================================================
// SurfaceLoads Tests
// ============================================================================

TEST(SurfaceLoadsTest, DefaultState) {
    SurfaceLoads sl;
    EXPECT_EQ(sl.surface_index, 0);
    EXPECT_TRUE(sl.name.empty());
    EXPECT_NEAR(sl.CL, 0.0, kTol);
    EXPECT_NEAR(sl.CD, 0.0, kTol);
    EXPECT_NEAR(sl.CY, 0.0, kTol);
    EXPECT_NEAR(sl.area, 0.0, kTol);
    EXPECT_TRUE(sl.strip_loads.empty());
}

}  // namespace
}  // namespace pegasus
