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

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace pegasus {

namespace {

// 1/(4*pi)
constexpr Scalar kPi4Inv = 0.079577471545947667884;

/// Velocity triple returned by the sweep kernel below.
struct KernelVelocity {
    Scalar u;
    Scalar v;
    Scalar w;
};

/// Branch-free horseshoe-vortex kernel for the assembly sweeps. Mirrors
/// ComputeHorseshoeVortexVelocity expression-for-expression (same operation
/// order, so results match to the last bit); the control-flow guards become
/// value selects so the field-point sweep can be inlined and vectorized.
/// Guarded terms may speculatively evaluate to inf/NaN when their guard is
/// false — the select discards them before they reach the accumulators.
inline KernelVelocity HorseshoeSweepKernel(Scalar fx, Scalar fy, Scalar fz,
                                           Scalar lx, Scalar ly, Scalar lz,
                                           Scalar rx, Scalar ry, Scalar rz,
                                           Scalar prandtl_glauert_factor,
                                           Scalar core_4, bool include_bound) {
    const Scalar ax = (lx - fx) / prandtl_glauert_factor;
    const Scalar ay = ly - fy;
    const Scalar az = lz - fz;
    const Scalar bx = (rx - fx) / prandtl_glauert_factor;
    const Scalar by = ry - fy;
    const Scalar bz = rz - fz;

    const Scalar asq = ax * ax + ay * ay + az * az;
    const Scalar bsq = bx * bx + by * by + bz * bz;
    const Scalar amag = std::sqrt(asq);
    const Scalar bmag = std::sqrt(bsq);

    Scalar U = 0.0;
    Scalar V = 0.0;
    Scalar W = 0.0;

    // Transverse bound leg (Leishman R^4 core model).
    const Scalar axb_x = ay * bz - az * by;
    const Scalar axb_y = az * bx - ax * bz;
    const Scalar axb_z = ax * by - ay * bx;
    const Scalar axb_sq = axb_x * axb_x + axb_y * axb_y + axb_z * axb_z;
    const Scalar adb = ax * bx + ay * by + az * bz;
    const Scalar alsq = asq + bsq - 2.0 * adb;
    const Scalar t_bound =
        ((bsq - adb) / std::sqrt(std::sqrt(bsq * bsq + core_4)) +
         (asq - adb) / std::sqrt(std::sqrt(asq * asq + core_4))) /
        std::sqrt(axb_sq * axb_sq + alsq * alsq * core_4);
    const bool bound_active =
        include_bound && (amag * bmag > kEpsilon) && (axb_sq > kEpsilon);
    const Scalar tb = bound_active ? t_bound : Scalar(0.0);
    U = axb_x * tb;
    V = axb_y * tb;
    W = axb_z * tb;

    // Trailing leg attached to A (semi-infinite along +x).
    const Scalar axisq = ay * ay + az * az;
    const Scalar t_a_raw =
        -(1.0 - ax / amag) / std::sqrt(axisq * axisq + core_4);
    const Scalar t_a = (amag > kEpsilon) ? t_a_raw : Scalar(0.0);
    V += az * t_a;
    W -= ay * t_a;

    // Trailing leg attached to B (semi-infinite along +x).
    const Scalar bxisq = by * by + bz * bz;
    const Scalar t_b_raw =
        (1.0 - bx / bmag) / std::sqrt(bxisq * bxisq + core_4);
    const Scalar t_b = (bmag > kEpsilon) ? t_b_raw : Scalar(0.0);
    V += bz * t_b;
    W -= by * t_b;

    return {U * kPi4Inv / prandtl_glauert_factor, V * kPi4Inv, W * kPi4Inv};
}

}  // namespace

Vec3 ComputeHorseshoeVortexVelocity(const Vec3& field_point,
                                    const Vec3& bound_left,
                                    const Vec3& bound_right,
                                    Scalar prandtl_glauert_factor,
                                    Scalar core_radius, bool include_bound) {
    // Prandtl-Glauert transformed coordinates
    // A = vector from bound_left to field_point
    // B = vector from bound_right to field_point
    Vec3 A, B;
    A.x() = (bound_left.x() - field_point.x()) / prandtl_glauert_factor;
    A.y() = bound_left.y() - field_point.y();
    A.z() = bound_left.z() - field_point.z();

    B.x() = (bound_right.x() - field_point.x()) / prandtl_glauert_factor;
    B.y() = bound_right.y() - field_point.y();
    B.z() = bound_right.z() - field_point.z();

    Scalar Asq = A.squaredNorm();
    Scalar Bsq = B.squaredNorm();
    Scalar Amag = std::sqrt(Asq);
    Scalar Bmag = std::sqrt(Bsq);

    Scalar core_sq = core_radius * core_radius;
    Scalar core_4 = core_sq * core_sq;

    Scalar U = 0.0, V = 0.0, W = 0.0;

    // Contribution from the transverse bound leg
    if (include_bound && Amag * Bmag > kEpsilon) {
        Vec3 AxB = A.cross(B);
        Scalar AxBsq = AxB.squaredNorm();

        if (AxBsq > kEpsilon) {
            Scalar AdB = A.dot(B);
            Scalar ALsq = Asq + Bsq - 2.0 * AdB;

            // Leishman R^4 core model (regularizes singularity)
            // Vtan = Gam/2*pi . r/sqrt(r^4 + rcore^4)
            Scalar T =
                ((Bsq - AdB) / std::sqrt(std::sqrt(Bsq * Bsq + core_4)) +
                 (Asq - AdB) / std::sqrt(std::sqrt(Asq * Asq + core_4))) /
                std::sqrt(AxBsq * AxBsq + ALsq * ALsq * core_4);

            U = AxB.x() * T;
            V = AxB.y() * T;
            W = AxB.z() * T;
        }
    }

    // Contribution from trailing leg attached to A (semi-infinite along +x)
    if (Amag > kEpsilon) {
        Scalar Axisq = A.y() * A.y() + A.z() * A.z();
        Scalar Adx = A.x();

        // Leishman core model for trailing leg
        Scalar T = -(1.0 - Adx / Amag) / std::sqrt(Axisq * Axisq + core_4);

        V += A.z() * T;
        W -= A.y() * T;
    }

    // Contribution from trailing leg attached to B (semi-infinite along +x)
    if (Bmag > kEpsilon) {
        Scalar Bxisq = B.y() * B.y() + B.z() * B.z();
        Scalar Bdx = B.x();

        // Leishman core model for trailing leg
        Scalar T = (1.0 - Bdx / Bmag) / std::sqrt(Bxisq * Bxisq + core_4);

        V += B.z() * T;
        W -= B.y() * T;
    }

    // Apply 1/(4*pi) factor and Prandtl-Glauert correction to U
    return {U * kPi4Inv / prandtl_glauert_factor, V * kPi4Inv, W * kPi4Inv};
}

InfluenceMatrixBuilder::InfluenceMatrixBuilder(const VortexLattice& mesh)
    : mesh_(mesh) {}

void InfluenceMatrixBuilder::SetMachNumber(Scalar mach) {
    mach_number_ = mach;
    if (mach < 1.0) {
        prandtl_glauert_factor_ = std::sqrt(1.0 - mach * mach);
    } else {
        prandtl_glauert_factor_ = 1.0;  // Invalid but avoid NaN
    }
}

Scalar InfluenceMatrixBuilder::ComputeCoreRadius(Index vortex_index,
                                                 bool same_component) const {
    const HorseshoeVortex& vor = mesh_.GetVortex(vortex_index);

    // Panel width (spanwise extent)
    Scalar panel_width =
        (vor.bound_vortex_right - vor.bound_vortex_left).norm();

    // Default small core based on panel width
    Scalar core = 0.0001 * panel_width;

    // Use larger core for different components to avoid near-field
    // singularities
    if (!same_component) {
        core = std::max(core_radius_chord_fraction_ * vor.local_chord,
                        core_radius_width_fraction_ * panel_width);
    }

    return core;
}

Vec3 InfluenceMatrixBuilder::ComputeSingleVortexVelocity(
    const Vec3& field_point, Index vortex_index, bool include_bound) const {
    const HorseshoeVortex& vor = mesh_.GetVortex(vortex_index);

    // For self-influence, use same component
    Scalar core = ComputeCoreRadius(vortex_index, true);

    return ComputeHorseshoeVortexVelocity(
        field_point, vor.bound_vortex_left, vor.bound_vortex_right,
        prandtl_glauert_factor_, core, include_bound);
}

Vec3 InfluenceMatrixBuilder::ComputeVelocityWithYImage(
    const Vec3& field_point, Index vortex_index,
    [[maybe_unused]] bool include_bound, Scalar core_radius) const {
    const HorseshoeVortex& vor = mesh_.GetVortex(vortex_index);

    Scalar y_offset = 2.0 * symmetry_.y_symmetry_plane;
    auto y_sign = static_cast<Scalar>(static_cast<int>(symmetry_.y_symmetry));

    // Y-image vortex: reflect about y-plane, swap left/right to maintain sign
    Vec3 image_left(vor.bound_vortex_right.x(),
                    y_offset - vor.bound_vortex_right.y(),
                    vor.bound_vortex_right.z());
    Vec3 image_right(vor.bound_vortex_left.x(),
                     y_offset - vor.bound_vortex_left.y(),
                     vor.bound_vortex_left.z());

    // Check if control point is at image vortex midpoint (skip bound for
    // symmetric case)
    bool include_image_bound = true;
    if (symmetry_.y_symmetry == SymmetryType::kSymmetric) {
        Vec3 image_mid = 0.5 * (image_left + image_right);
        if ((field_point - image_mid).squaredNorm() < kEpsilon) {
            include_image_bound = false;
        }
    }

    Vec3 vel = ComputeHorseshoeVortexVelocity(
        field_point, image_left, image_right, prandtl_glauert_factor_,
        core_radius, include_image_bound);

    return y_sign * vel;
}

Vec3 InfluenceMatrixBuilder::ComputeVelocityWithZImage(
    const Vec3& field_point, Index vortex_index,
    [[maybe_unused]] bool include_bound, Scalar core_radius) const {
    const HorseshoeVortex& vor = mesh_.GetVortex(vortex_index);

    Scalar z_offset = 2.0 * symmetry_.z_symmetry_plane;
    auto z_sign = static_cast<Scalar>(static_cast<int>(symmetry_.z_symmetry));

    // Z-image vortex: reflect about z-plane, swap left/right
    Vec3 image_left(vor.bound_vortex_right.x(), vor.bound_vortex_right.y(),
                    z_offset - vor.bound_vortex_right.z());
    Vec3 image_right(vor.bound_vortex_left.x(), vor.bound_vortex_left.y(),
                     z_offset - vor.bound_vortex_left.z());

    Vec3 vel = ComputeHorseshoeVortexVelocity(
        field_point, image_left, image_right, prandtl_glauert_factor_,
        core_radius, true);

    return z_sign * vel;
}

MatX InfluenceMatrixBuilder::BuildVelocityInfluenceMatrix() const {
    const Index n = mesh_.GetNumVortices();
    MatX wc_gam(3 * n, n);

    // Contiguous per-field-point data for the inner sweeps.
    std::vector<Scalar> fx(static_cast<std::size_t>(n));
    std::vector<Scalar> fy(static_cast<std::size_t>(n));
    std::vector<Scalar> fz(static_cast<std::size_t>(n));
    std::vector<Index> component(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh_.GetVortex(i);
        const auto idx = static_cast<std::size_t>(i);
        fx[idx] = vor.control_point.x();
        fy[idx] = vor.control_point.y();
        fz[idx] = vor.control_point.z();
        component[idx] = vor.component_index;
    }

    const bool has_y = symmetry_.y_symmetry != SymmetryType::kNone;
    const bool has_z = symmetry_.z_symmetry != SymmetryType::kNone;
    const bool y_symmetric = symmetry_.y_symmetry == SymmetryType::kSymmetric;
    const Scalar y_offset = 2.0 * symmetry_.y_symmetry_plane;
    const Scalar z_offset = 2.0 * symmetry_.z_symmetry_plane;
    const auto y_sign =
        static_cast<Scalar>(static_cast<int>(symmetry_.y_symmetry));
    const auto z_sign =
        static_cast<Scalar>(static_cast<int>(symmetry_.z_symmetry));
    const Scalar yz_sign = y_sign * z_sign;
    const Scalar pg = prandtl_glauert_factor_;

    // Source-major sweep: all per-source data (endpoints, image endpoints,
    // core radii) is hoisted out of the field-point loop, and each source
    // writes one contiguous column of the column-major matrix. Columns are
    // independent, so the outer loop parallelizes deterministically.
#ifdef PEGASUS_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (Index j = 0; j < n; ++j) {
        const HorseshoeVortex& vor_j = mesh_.GetVortex(j);
        const Vec3 left = vor_j.bound_vortex_left;
        const Vec3 right = vor_j.bound_vortex_right;

        // Per-source core radii for same/different components (previously
        // recomputed per (i, j) pair inside the inner loop).
        const Scalar panel_width = (right - left).norm();
        const Scalar core_same = 0.0001 * panel_width;
        const Scalar core_diff =
            std::max(core_radius_chord_fraction_ * vor_j.local_chord,
                     core_radius_width_fraction_ * panel_width);
        const Scalar core_same_sq = core_same * core_same;
        const Scalar core4_same = core_same_sq * core_same_sq;
        const Scalar core_diff_sq = core_diff * core_diff;
        const Scalar core4_diff = core_diff_sq * core_diff_sq;
        const Index comp_j = vor_j.component_index;

        // Y-image vortex: reflect about y-plane, swap left/right.
        const Vec3 yl(right.x(), y_offset - right.y(), right.z());
        const Vec3 yr(left.x(), y_offset - left.y(), left.z());
        const Vec3 y_mid = 0.5 * (yl + yr);
        // Z-image vortex: reflect about z-plane, swap left/right.
        const Vec3 zl(right.x(), right.y(), z_offset - right.z());
        const Vec3 zr(left.x(), left.y(), z_offset - left.z());
        // Y-Z double image (no swap, matching the original assembly).
        const Vec3 yzl(left.x(), y_offset - left.y(), z_offset - left.z());
        const Vec3 yzr(right.x(), y_offset - right.y(), z_offset - right.z());

        Scalar* col = wc_gam.data() + 3 * n * j;
        // The kernel is branch-free (guards are value selects), so this
        // sweep is legal to auto-vectorize; whether that pays is up to the
        // target's cost model (on Apple Silicon the scalar sqrt/div
        // pipelines win, so it is deliberately not forced).
        for (Index i = 0; i < n; ++i) {
            const auto idx = static_cast<std::size_t>(i);
            const Scalar core_4 =
                (component[idx] == comp_j) ? core4_same : core4_diff;

            // Real vortex contribution (skip bound for self-influence).
            KernelVelocity vel = HorseshoeSweepKernel(
                fx[idx], fy[idx], fz[idx], left.x(), left.y(), left.z(),
                right.x(), right.y(), right.z(), pg, core_4, i != j);

            if (has_y) {
                // Skip the image bound leg if the control point sits at the
                // image vortex midpoint (symmetric case only).
                bool include_image_bound = true;
                if (y_symmetric) {
                    const Scalar dx = fx[idx] - y_mid.x();
                    const Scalar dy = fy[idx] - y_mid.y();
                    const Scalar dz = fz[idx] - y_mid.z();
                    include_image_bound =
                        (dx * dx + dy * dy + dz * dz) >= kEpsilon;
                }
                const KernelVelocity vy = HorseshoeSweepKernel(
                    fx[idx], fy[idx], fz[idx], yl.x(), yl.y(), yl.z(), yr.x(),
                    yr.y(), yr.z(), pg, core_4, include_image_bound);
                vel.u += y_sign * vy.u;
                vel.v += y_sign * vy.v;
                vel.w += y_sign * vy.w;
            }

            if (has_z) {
                const KernelVelocity vz = HorseshoeSweepKernel(
                    fx[idx], fy[idx], fz[idx], zl.x(), zl.y(), zl.z(), zr.x(),
                    zr.y(), zr.z(), pg, core_4, true);
                vel.u += z_sign * vz.u;
                vel.v += z_sign * vz.v;
                vel.w += z_sign * vz.w;

                if (has_y) {
                    const KernelVelocity vyz = HorseshoeSweepKernel(
                        fx[idx], fy[idx], fz[idx], yzl.x(), yzl.y(), yzl.z(),
                        yzr.x(), yzr.y(), yzr.z(), pg, core_4, true);
                    vel.u += yz_sign * vyz.u;
                    vel.v += yz_sign * vyz.v;
                    vel.w += yz_sign * vyz.w;
                }
            }

            col[3 * i + 0] = vel.u;
            col[3 * i + 1] = vel.v;
            col[3 * i + 2] = vel.w;
        }
    }

    return wc_gam;
}

MatX InfluenceMatrixBuilder::BuildNormalwashMatrix() const {
    return BuildNormalwashMatrix(nullptr);
}

MatX InfluenceMatrixBuilder::BuildNormalwashMatrix(
    MatX* velocity_matrix) const {
    MatX wc_gam = BuildVelocityInfluenceMatrix();

    const Index n = mesh_.GetNumVortices();
    MatX aic(n, n);

    // Contiguous per-control-point normals for the column sweeps.
    std::vector<Scalar> nx(static_cast<std::size_t>(n));
    std::vector<Scalar> ny(static_cast<std::size_t>(n));
    std::vector<Scalar> nz(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i) {
        const Vec3& normal = mesh_.GetVortex(i).normal_at_control_point;
        const auto idx = static_cast<std::size_t>(i);
        nx[idx] = normal.x();
        ny[idx] = normal.y();
        nz[idx] = normal.z();
    }

    // Project velocity influence onto panel normals, column by column so
    // both matrices are traversed contiguously (column-major storage).
    for (Index j = 0; j < n; ++j) {
        const Scalar* col = wc_gam.data() + 3 * n * j;
        Scalar* aic_col = aic.data() + n * j;
        for (Index i = 0; i < n; ++i) {
            const auto idx = static_cast<std::size_t>(i);
            aic_col[i] = col[3 * i + 0] * nx[idx] + col[3 * i + 1] * ny[idx] +
                         col[3 * i + 2] * nz[idx];
        }
    }

    if (velocity_matrix != nullptr) {
        *velocity_matrix = std::move(wc_gam);
    }

    return aic;
}

Vec3 InfluenceMatrixBuilder::ComputeInducedVelocity(
    const Vec3& field_point, const VecX& circulation) const {
    Vec3 velocity = Vec3::Zero();
    Index n = mesh_.GetNumVortices();

    for (Index j = 0; j < n; ++j) {
        const HorseshoeVortex& vor_j = mesh_.GetVortex(j);
        Scalar core = ComputeCoreRadius(j, false);

        Vec3 vel = ComputeHorseshoeVortexVelocity(
            field_point, vor_j.bound_vortex_left, vor_j.bound_vortex_right,
            prandtl_glauert_factor_, core, true);

        // Add symmetry images
        if (symmetry_.y_symmetry != SymmetryType::kNone) {
            vel += ComputeVelocityWithYImage(field_point, j, true, core);
        }
        if (symmetry_.z_symmetry != SymmetryType::kNone) {
            vel += ComputeVelocityWithZImage(field_point, j, true, core);
        }

        velocity += circulation(j) * vel;
    }

    return velocity;
}

}  // namespace pegasus
