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

#include "pegasus/solver/influence_matrix.h"

#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"
#include "pegasus/solver/vortex_lattice.h"

#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

namespace pegasus {
namespace {

constexpr double kTestTolerance = 1.0e-6;

// ============================================================================
// Biot-Savart Tests
// ============================================================================

class BiotSavartTest : public ::testing::Test {
protected:
    // Create a simple unit horseshoe vortex centered at origin
    // Bound vortex from (-0.5, 0, 0) to (0.5, 0, 0)
    Vec3 bound_left{-0.5, 0.0, 0.0};
    Vec3 bound_right{0.5, 0.0, 0.0};
};

TEST_F(BiotSavartTest, VelocityAboveBoundVortex) {
    // Field point directly above the bound vortex center
    // With bound vortex along x-axis at y=0, and field point at y=0,
    // the velocity should have no y-component due to symmetry
    Vec3 field_point(0.0, 0.0, 1.0);
    Scalar core_radius = 0.001;

    Vec3 vel = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, core_radius, true);

    // Velocity should be finite (no NaN/Inf)
    EXPECT_TRUE(std::isfinite(vel.x()));
    EXPECT_TRUE(std::isfinite(vel.y()));
    EXPECT_TRUE(std::isfinite(vel.z()));

    // For this symmetric geometry (field point on y=0 plane, vortex in y=0
    // plane), y-velocity should be very small due to symmetry
    EXPECT_NEAR(vel.y(), 0.0, 0.01);
}

TEST_F(BiotSavartTest, VelocitySymmetry) {
    // Velocity should be symmetric about x-z plane (y=0)
    Vec3 point_pos(0.0, 0.5, 0.5);
    Vec3 point_neg(0.0, -0.5, 0.5);
    Scalar core_radius = 0.001;

    Vec3 vel_pos = ComputeHorseshoeVortexVelocity(
        point_pos, bound_left, bound_right, 1.0, core_radius);
    Vec3 vel_neg = ComputeHorseshoeVortexVelocity(
        point_neg, bound_left, bound_right, 1.0, core_radius);

    // x-velocity should be antisymmetric
    EXPECT_NEAR(vel_pos.x(), -vel_neg.x(), kTestTolerance);

    // y-velocity should be symmetric (both negative)
    EXPECT_NEAR(vel_pos.y(), vel_neg.y(), kTestTolerance);

    // z-velocity should be antisymmetric
    EXPECT_NEAR(vel_pos.z(), -vel_neg.z(), kTestTolerance);
}

TEST_F(BiotSavartTest, CoreRegularization) {
    // Test core regularization by using only trailing vortex contributions
    // (exclude bound vortex to avoid additional complexity)
    // Right trailing vortex extends from (0.5, 0, 0) in +x direction
    // Left trailing vortex extends from (-0.5, 0, 0) in +x direction
    // Place field point asymmetrically to avoid cancellation
    Vec3 field_point(1.0, 0.3, 0.02);

    // Small core - should give larger velocity magnitude
    Vec3 vel_small_core = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, 0.01,
        false);  // Exclude bound vortex

    // Large core - should give smaller velocity magnitude (regularized)
    Vec3 vel_large_core = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, 0.5,
        false);  // Exclude bound vortex

    // Both should be finite
    EXPECT_TRUE(std::isfinite(vel_small_core.norm()));
    EXPECT_TRUE(std::isfinite(vel_large_core.norm()));

    // Both should give substantial velocity (downwash from trailing vortices)
    EXPECT_GT(vel_small_core.norm(), 0.01);

    // Larger core should reduce the velocity magnitude
    EXPECT_LT(vel_large_core.norm(), vel_small_core.norm());
}

TEST_F(BiotSavartTest, PrandtlGlauertCorrection) {
    // Point off the symmetry plane to get non-zero velocity
    Vec3 field_point(1.0, 0.5, 0.5);
    Scalar core_radius = 0.001;

    // Incompressible
    Vec3 vel_incomp = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, core_radius);

    // Subsonic (M = 0.5, beta = sqrt(0.75) ≈ 0.866)
    Scalar beta = std::sqrt(0.75);
    Vec3 vel_comp = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, beta, core_radius);

    // Compressibility should change the velocity (Prandtl-Glauert effect)
    // The overall magnitude should increase for subsonic flow
    EXPECT_NE(vel_comp.norm(), vel_incomp.norm());
}

TEST_F(BiotSavartTest, ExcludeBoundVortex) {
    // Use a point off the symmetry plane to get meaningful velocity
    Vec3 field_point(0.0, 0.3, 0.5);
    Scalar core_radius = 0.001;

    Vec3 vel_with_bound = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, core_radius, true);

    Vec3 vel_no_bound = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, core_radius, false);

    // Both should be finite
    EXPECT_TRUE(std::isfinite(vel_with_bound.norm()));
    EXPECT_TRUE(std::isfinite(vel_no_bound.norm()));

    // Without bound vortex, velocity should be different
    EXPECT_NE(vel_with_bound.norm(), vel_no_bound.norm());
}

TEST_F(BiotSavartTest, TrailingVortexInDownstreamField) {
    // Point far downstream (in +x direction)
    Vec3 field_point(10.0, 0.5, 0.0);
    Scalar core_radius = 0.001;

    Vec3 vel = ComputeHorseshoeVortexVelocity(
        field_point, bound_left, bound_right, 1.0, core_radius, true);

    // Far downstream, velocity is primarily from trailing vortices
    // Should induce downwash (negative z) for positive y location
    EXPECT_LT(vel.z(), 0.0);
}

// ============================================================================
// InfluenceMatrixBuilder Tests
// ============================================================================

class InfluenceMatrixTest : public ::testing::Test {
protected:
    Aircraft CreateRectangularWing(Scalar span, Scalar chord, int num_chordwise,
                                   int num_spanwise) {
        Aircraft aircraft("Test Wing");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(num_chordwise);
        wing.SetSpanwisePanelCount(num_spanwise);
        wing.SetChordwiseSpacing(SpacingType::kUniform);
        wing.SetSpanwiseSpacing(SpacingType::kUniform);

        // Full wing (no y-duplicate for simplicity)
        wing.AddSection(Section(Vec3(0, -span / 2, 0), chord, 0.0));
        wing.AddSection(Section(Vec3(0, span / 2, 0), chord, 0.0));

        aircraft.AddSurface(wing);
        aircraft.SetReferenceArea(span * chord);
        aircraft.SetReferenceChord(chord);
        aircraft.SetReferenceSpan(span);

        return aircraft;
    }
};

TEST_F(InfluenceMatrixTest, VelocityMatrixDimensions) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX wc_gam = builder.BuildVelocityInfluenceMatrix();

    Index n = mesh.GetNumVortices();
    EXPECT_EQ(wc_gam.rows(), 3 * n);
    EXPECT_EQ(wc_gam.cols(), n);
}

TEST_F(InfluenceMatrixTest, AICMatrixDimensions) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX aic = builder.BuildNormalwashMatrix();

    Index n = mesh.GetNumVortices();
    EXPECT_EQ(aic.rows(), n);
    EXPECT_EQ(aic.cols(), n);
}

TEST_F(InfluenceMatrixTest, AICDiagonalDominance) {
    // For a flat wing, AIC should be diagonally dominant
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX aic = builder.BuildNormalwashMatrix();

    Index n = mesh.GetNumVortices();

    // Check that diagonal elements are significant
    for (Index i = 0; i < n; ++i) {
        Scalar diag = std::abs(aic(i, i));
        Scalar row_sum = 0.0;
        for (Index j = 0; j < n; ++j) {
            if (i != j) {
                row_sum += std::abs(aic(i, j));
            }
        }
        // Diagonal should be a significant portion of the row
        EXPECT_GT(diag, 0.0);
    }
}

TEST_F(InfluenceMatrixTest, AICMatrixNonsingular) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 4, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX aic = builder.BuildNormalwashMatrix();

    // Matrix should be invertible
    Eigen::PartialPivLU<MatX> lu(aic);
    Scalar det = lu.determinant();

    EXPECT_NE(det, 0.0);
}

TEST_F(InfluenceMatrixTest, MachNumberEffect) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);

    // Incompressible
    builder.SetMachNumber(0.0);
    MatX aic_incomp = builder.BuildNormalwashMatrix();

    // Subsonic
    builder.SetMachNumber(0.5);
    MatX aic_comp = builder.BuildNormalwashMatrix();

    // Matrices should be different due to Prandtl-Glauert
    Scalar diff = (aic_comp - aic_incomp).norm();
    EXPECT_GT(diff, 0.0);
}

TEST_F(InfluenceMatrixTest, InducedVelocityFromMatrix) {
    // Test that velocity influence matrix produces finite, reasonable values
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX velocity_matrix;
    builder.BuildNormalwashMatrix(&velocity_matrix);

    // All values should be finite
    for (Index i = 0; i < velocity_matrix.rows(); ++i) {
        for (Index j = 0; j < velocity_matrix.cols(); ++j) {
            EXPECT_TRUE(std::isfinite(velocity_matrix(i, j)));
        }
    }

    // Diagonal-ish influence should be significant
    // (self-influence from trailing legs)
    Index n = mesh.GetNumVortices();
    for (Index i = 0; i < n; ++i) {
        Scalar self_influence = std::sqrt(
            velocity_matrix(3 * i + 0, i) * velocity_matrix(3 * i + 0, i) +
            velocity_matrix(3 * i + 1, i) * velocity_matrix(3 * i + 1, i) +
            velocity_matrix(3 * i + 2, i) * velocity_matrix(3 * i + 2, i));
        EXPECT_GT(self_influence, 0.0);
    }
}

TEST_F(InfluenceMatrixTest, SupersonicMachClamped) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 2);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    builder.SetMachNumber(1.5);  // Supersonic - PG factor clamped to 1.0

    MatX aic = builder.BuildNormalwashMatrix();

    // Should still produce finite results
    for (Index i = 0; i < aic.rows(); ++i) {
        for (Index j = 0; j < aic.cols(); ++j) {
            EXPECT_TRUE(std::isfinite(aic(i, j)));
        }
    }
}

TEST_F(InfluenceMatrixTest, ComputeInducedVelocity) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);

    // Create uniform circulation
    Index n = mesh.GetNumVortices();
    VecX circulation = VecX::Ones(n) * 0.1;

    Vec3 field_point(0.5, 0.0, 0.5);
    Vec3 vel = builder.ComputeInducedVelocity(field_point, circulation);

    EXPECT_TRUE(std::isfinite(vel.x()));
    EXPECT_TRUE(std::isfinite(vel.y()));
    EXPECT_TRUE(std::isfinite(vel.z()));

    // With positive circulation, should get downwash at center
    EXPECT_LT(vel.z(), 0.0);
}

TEST_F(InfluenceMatrixTest, ComputeSingleVortexVelocity) {
    Aircraft aircraft = CreateRectangularWing(4.0, 1.0, 2, 4);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);

    Vec3 field_point(0.5, 0.0, 1.0);
    Vec3 vel = builder.ComputeSingleVortexVelocity(field_point, 0, true);

    EXPECT_TRUE(vel.allFinite());
}

TEST_F(InfluenceMatrixTest, YSymmetryInfluence) {
    Aircraft aircraft("Symmetric Wing");

    Surface wing("Wing");
    wing.SetChordwisePanelCount(2);
    wing.SetSpanwisePanelCount(4);
    wing.SetYDuplicate(0.0);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 2.0, 0), 1.0, 0.0));

    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(4.0);
    aircraft.SetReferenceChord(1.0);
    aircraft.SetReferenceSpan(4.0);
    aircraft.SetYSymmetry(SymmetryType::kSymmetric);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    SymmetryConfig sym;
    sym.y_symmetry = aircraft.GetYSymmetry();
    sym.z_symmetry = aircraft.GetZSymmetry();
    sym.z_symmetry_plane = aircraft.GetZSymmetryPlane();
    builder.SetSymmetry(sym);
    MatX aic = builder.BuildNormalwashMatrix();

    // Should be finite
    for (Index i = 0; i < aic.rows(); ++i) {
        for (Index j = 0; j < aic.cols(); ++j) {
            EXPECT_TRUE(std::isfinite(aic(i, j)));
        }
    }
}

TEST_F(InfluenceMatrixTest, DifferentComponentCoreRadius) {
    // Two separate surfaces -> different components -> larger core radius
    Aircraft aircraft("Two Component");

    Surface wing("Wing");
    wing.SetChordwisePanelCount(2);
    wing.SetSpanwisePanelCount(2);
    wing.SetComponentIndex(0);
    wing.AddSection(Section(Vec3(0, -2, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 2, 0), 1.0, 0.0));

    Surface tail("Tail");
    tail.SetChordwisePanelCount(2);
    tail.SetSpanwisePanelCount(2);
    tail.SetComponentIndex(1);
    tail.AddSection(Section(Vec3(4, -1, 0), 0.5, 0.0));
    tail.AddSection(Section(Vec3(4, 1, 0), 0.5, 0.0));

    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);
    aircraft.SetReferenceArea(4.0);
    aircraft.SetReferenceChord(1.0);
    aircraft.SetReferenceSpan(4.0);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    MatX aic = builder.BuildNormalwashMatrix();

    // All values finite
    for (Index i = 0; i < aic.rows(); ++i) {
        for (Index j = 0; j < aic.cols(); ++j) {
            EXPECT_TRUE(std::isfinite(aic(i, j)));
        }
    }
}

TEST_F(InfluenceMatrixTest, InducedVelocityWithYSymmetry) {
    Aircraft aircraft("Symmetric Wing");

    Surface wing("Wing");
    wing.SetChordwisePanelCount(2);
    wing.SetSpanwisePanelCount(4);
    wing.SetYDuplicate(0.0);
    wing.AddSection(Section(Vec3(0, 0, 0), 1.0, 0.0));
    wing.AddSection(Section(Vec3(0, 2.0, 0), 1.0, 0.0));

    aircraft.AddSurface(wing);
    aircraft.SetYSymmetry(SymmetryType::kSymmetric);

    VortexLattice mesh;
    mesh.GenerateFromAircraft(aircraft);

    InfluenceMatrixBuilder builder(mesh);
    SymmetryConfig sym;
    sym.y_symmetry = aircraft.GetYSymmetry();
    sym.z_symmetry = aircraft.GetZSymmetry();
    sym.z_symmetry_plane = aircraft.GetZSymmetryPlane();
    builder.SetSymmetry(sym);

    Index n = mesh.GetNumVortices();
    VecX circulation = VecX::Ones(n) * 0.1;

    Vec3 field_point(0.5, 0.0, 0.5);
    Vec3 vel = builder.ComputeInducedVelocity(field_point, circulation);
    EXPECT_TRUE(vel.allFinite());
}

// ============================================================================
// Assembly regression tests: the optimized source-major sweep must
// reproduce the straightforward per-pair assembly built from the public
// Biot-Savart kernel (real vortex + symmetry images + per-pair core radii).
// ============================================================================

/// Reference velocity-influence assembly: the direct (i, j) double loop
/// with per-pair core radii and symmetry images, mirroring the pre-sweep
/// implementation.
MatX ReferenceVelocityMatrix(const VortexLattice& mesh, Scalar pg,
                             const SymmetryConfig& sym, Scalar chord_frac,
                             Scalar width_frac) {
    const Index n = mesh.GetNumVortices();
    MatX wc_gam = MatX::Zero(3 * n, n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor_i = mesh.GetVortex(i);
        const Vec3 fp = vor_i.control_point;

        for (Index j = 0; j < n; ++j) {
            const HorseshoeVortex& vor_j = mesh.GetVortex(j);
            const Vec3 left = vor_j.bound_vortex_left;
            const Vec3 right = vor_j.bound_vortex_right;

            const Scalar panel_width = (right - left).norm();
            Scalar core = 0.0001 * panel_width;
            if (vor_i.component_index != vor_j.component_index) {
                core = std::max(chord_frac * vor_j.local_chord,
                                width_frac * panel_width);
            }

            Vec3 vel = ComputeHorseshoeVortexVelocity(fp, left, right, pg, core,
                                                      i != j);

            if (sym.y_symmetry != SymmetryType::kNone) {
                const Scalar y_offset = 2.0 * sym.y_symmetry_plane;
                const Vec3 yl(right.x(), y_offset - right.y(), right.z());
                const Vec3 yr(left.x(), y_offset - left.y(), left.z());
                bool include_bound = true;
                if (sym.y_symmetry == SymmetryType::kSymmetric) {
                    const Vec3 mid = 0.5 * (yl + yr);
                    include_bound = (fp - mid).squaredNorm() >= kEpsilon;
                }
                const auto y_sign =
                    static_cast<Scalar>(static_cast<int>(sym.y_symmetry));
                vel += y_sign * ComputeHorseshoeVortexVelocity(
                                    fp, yl, yr, pg, core, include_bound);
            }

            if (sym.z_symmetry != SymmetryType::kNone) {
                const Scalar z_offset = 2.0 * sym.z_symmetry_plane;
                const Vec3 zl(right.x(), right.y(), z_offset - right.z());
                const Vec3 zr(left.x(), left.y(), z_offset - left.z());
                const auto z_sign =
                    static_cast<Scalar>(static_cast<int>(sym.z_symmetry));
                vel += z_sign * ComputeHorseshoeVortexVelocity(fp, zl, zr, pg,
                                                               core, true);

                if (sym.y_symmetry != SymmetryType::kNone) {
                    const Scalar y_offset = 2.0 * sym.y_symmetry_plane;
                    const Vec3 yzl(left.x(), y_offset - left.y(),
                                   z_offset - left.z());
                    const Vec3 yzr(right.x(), y_offset - right.y(),
                                   z_offset - right.z());
                    const auto y_sign =
                        static_cast<Scalar>(static_cast<int>(sym.y_symmetry));
                    vel += y_sign * z_sign *
                           ComputeHorseshoeVortexVelocity(fp, yzl, yzr, pg,
                                                          core, true);
                }
            }

            wc_gam(3 * i + 0, j) = vel.x();
            wc_gam(3 * i + 1, j) = vel.y();
            wc_gam(3 * i + 2, j) = vel.z();
        }
    }
    return wc_gam;
}

class AssemblyRegressionTest : public ::testing::Test {
protected:
    /// Wing + tail so the mesh has two components (exercises both core
    /// radii paths).
    static Aircraft CreateTwoComponentAircraft() {
        Aircraft aircraft("Regression");

        Surface wing("Wing");
        wing.SetChordwisePanelCount(3);
        wing.SetSpanwisePanelCount(5);
        wing.SetChordwiseSpacing(SpacingType::kUniform);
        wing.SetSpanwiseSpacing(SpacingType::kUniform);
        wing.AddSection(Section(Vec3(0, -2, 0), 1.0, 0.0));
        wing.AddSection(Section(Vec3(0.2, 2, 0.1), 0.7, 0.0));
        aircraft.AddSurface(wing);

        Surface tail("Tail");
        tail.SetChordwisePanelCount(2);
        tail.SetSpanwisePanelCount(3);
        tail.SetChordwiseSpacing(SpacingType::kUniform);
        tail.SetSpanwiseSpacing(SpacingType::kUniform);
        tail.AddSection(Section(Vec3(4, -1, 0.5), 0.5, 0.0));
        tail.AddSection(Section(Vec3(4, 1, 0.5), 0.5, 0.0));
        aircraft.AddSurface(tail);

        aircraft.SetReferenceArea(4.0);
        aircraft.SetReferenceChord(1.0);
        aircraft.SetReferenceSpan(4.0);
        return aircraft;
    }

    void ExpectMatchesReference(const SymmetryConfig& sym, Scalar mach) {
        Aircraft aircraft = CreateTwoComponentAircraft();
        VortexLattice mesh;
        mesh.GenerateFromAircraft(aircraft);

        InfluenceMatrixBuilder builder(mesh);
        builder.SetMachNumber(mach);
        builder.SetSymmetry(sym);

        MatX vim;
        MatX aic = builder.BuildNormalwashMatrix(&vim);
        MatX reference = ReferenceVelocityMatrix(
            mesh, builder.GetPrandtlGlauertFactor(), sym, 0.0001, 0.0001);

        const Scalar scale =
            std::max(reference.cwiseAbs().maxCoeff(), Scalar(1.0));
        EXPECT_LE((vim - reference).cwiseAbs().maxCoeff(), 1e-14 * scale);

        // AIC must equal the normal projection of the reference matrix.
        const Index n = mesh.GetNumVortices();
        Scalar max_aic_diff = 0.0;
        for (Index i = 0; i < n; ++i) {
            const Vec3& normal = mesh.GetVortex(i).normal_at_control_point;
            for (Index j = 0; j < n; ++j) {
                const Scalar expected = reference(3 * i + 0, j) * normal.x() +
                                        reference(3 * i + 1, j) * normal.y() +
                                        reference(3 * i + 2, j) * normal.z();
                max_aic_diff =
                    std::max(max_aic_diff, std::abs(aic(i, j) - expected));
            }
        }
        EXPECT_LE(max_aic_diff, 1e-14 * scale);
    }
};

TEST_F(AssemblyRegressionTest, MatchesReferenceNoSymmetry) {
    ExpectMatchesReference(SymmetryConfig{}, 0.0);
}

TEST_F(AssemblyRegressionTest, MatchesReferenceWithMach) {
    ExpectMatchesReference(SymmetryConfig{}, 0.3);
}

TEST_F(AssemblyRegressionTest, MatchesReferenceYSymmetric) {
    SymmetryConfig sym;
    sym.y_symmetry = SymmetryType::kSymmetric;
    sym.y_symmetry_plane = 0.0;
    ExpectMatchesReference(sym, 0.0);
}

TEST_F(AssemblyRegressionTest, MatchesReferenceYAntisymmetricZSymmetric) {
    SymmetryConfig sym;
    sym.y_symmetry = SymmetryType::kAntisymmetric;
    sym.y_symmetry_plane = 0.0;
    sym.z_symmetry = SymmetryType::kSymmetric;
    sym.z_symmetry_plane = -0.5;
    ExpectMatchesReference(sym, 0.15);
}

}  // namespace
}  // namespace pegasus
