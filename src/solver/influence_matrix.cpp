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

#include "vlm/solver/influence_matrix.h"

#include <algorithm>
#include <cmath>

namespace vlm {

namespace {

// 1/(4*pi)
constexpr Scalar kPi4Inv = 0.079577471545947667884;

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
    return Vec3(U * kPi4Inv / prandtl_glauert_factor, V * kPi4Inv, W * kPi4Inv);
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
    Scalar y_sign = static_cast<Scalar>(static_cast<int>(symmetry_.y_symmetry));

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
    Scalar z_sign = static_cast<Scalar>(static_cast<int>(symmetry_.z_symmetry));

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
    Index n = mesh_.GetNumVortices();
    MatX wc_gam = MatX::Zero(3 * n, n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor_i = mesh_.GetVortex(i);
        Vec3 field_point = vor_i.control_point;

        for (Index j = 0; j < n; ++j) {
            const HorseshoeVortex& vor_j = mesh_.GetVortex(j);

            // Determine if on same component
            bool same_component =
                (vor_i.component_index == vor_j.component_index);
            Scalar core = ComputeCoreRadius(j, same_component);

            // Skip bound vortex contribution for self-influence
            bool include_bound = !(i == j);

            // Real vortex contribution
            Vec3 vel = ComputeHorseshoeVortexVelocity(
                field_point, vor_j.bound_vortex_left, vor_j.bound_vortex_right,
                prandtl_glauert_factor_, core, include_bound);

            // Y-symmetry image contribution
            if (symmetry_.y_symmetry != SymmetryType::kNone) {
                vel += ComputeVelocityWithYImage(field_point, j, true, core);
            }

            // Z-symmetry image contribution
            if (symmetry_.z_symmetry != SymmetryType::kNone) {
                vel += ComputeVelocityWithZImage(field_point, j, true, core);

                // Y-Z double image if both symmetries present
                if (symmetry_.y_symmetry != SymmetryType::kNone) {
                    Scalar y_offset = 2.0 * symmetry_.y_symmetry_plane;
                    Scalar z_offset = 2.0 * symmetry_.z_symmetry_plane;
                    Scalar yz_sign =
                        static_cast<Scalar>(
                            static_cast<int>(symmetry_.y_symmetry)) *
                        static_cast<Scalar>(
                            static_cast<int>(symmetry_.z_symmetry));

                    Vec3 image_left(vor_j.bound_vortex_left.x(),
                                    y_offset - vor_j.bound_vortex_left.y(),
                                    z_offset - vor_j.bound_vortex_left.z());
                    Vec3 image_right(vor_j.bound_vortex_right.x(),
                                     y_offset - vor_j.bound_vortex_right.y(),
                                     z_offset - vor_j.bound_vortex_right.z());

                    Vec3 yz_vel = ComputeHorseshoeVortexVelocity(
                        field_point, image_left, image_right,
                        prandtl_glauert_factor_, core, true);

                    vel += yz_sign * yz_vel;
                }
            }

            wc_gam(3 * i + 0, j) = vel.x();
            wc_gam(3 * i + 1, j) = vel.y();
            wc_gam(3 * i + 2, j) = vel.z();
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

    if (velocity_matrix != nullptr) {
        *velocity_matrix = wc_gam;
    }

    Index n = mesh_.GetNumVortices();
    MatX aic = MatX::Zero(n, n);

    // Project velocity influence onto panel normals
    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor_i = mesh_.GetVortex(i);
        Vec3 normal = vor_i.normal_at_control_point;

        for (Index j = 0; j < n; ++j) {
            aic(i, j) = wc_gam(3 * i + 0, j) * normal.x() +
                        wc_gam(3 * i + 1, j) * normal.y() +
                        wc_gam(3 * i + 2, j) * normal.z();
        }
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

}  // namespace vlm
