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

#include "pegasus/gradients/forward_jvp.h"

#include "pegasus/gradients/analytical_aic.h"
#include "pegasus/gradients/analytical_biot_savart.h"
#include "pegasus/gradients/analytical_forces.h"
#include "pegasus/gradients/mesh_sensitivity.h"

#include <cmath>
#include <vector>

namespace pegasus::gradients {

namespace {

/// Forward mode through normal = normalize(d1 x d2) where d1 = n11-n00,
/// d2 = n10-n01. Given perturbations of the four corner nodes, computes
/// the perturbation of the unit normal vector.
void ForwardNormalPerturbation(const double* n00, const double* n10,
                               const double* n01, const double* n11,
                               const double* dn00, const double* dn10,
                               const double* dn01, const double* dn11,
                               double* d_normal) {
    // Forward: d1 = n11 - n00, d2 = n10 - n01
    double d1[3] = {n11[0] - n00[0], n11[1] - n00[1], n11[2] - n00[2]};
    double d2[3] = {n10[0] - n01[0], n10[1] - n01[1], n10[2] - n01[2]};

    // raw = d1 x d2
    double raw[3];
    raw[0] = d1[1] * d2[2] - d1[2] * d2[1];
    raw[1] = d1[2] * d2[0] - d1[0] * d2[2];
    raw[2] = d1[0] * d2[1] - d1[1] * d2[0];

    double mag = std::sqrt(raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]);
    if (mag < 1e-12) {
        d_normal[0] = d_normal[1] = d_normal[2] = 0.0;
        return;
    }

    double inv_mag = 1.0 / mag;
    double normal[3] = {raw[0] * inv_mag, raw[1] * inv_mag, raw[2] * inv_mag};

    // Perturbations: dd1 = dn11 - dn00, dd2 = dn10 - dn01
    double dd1[3] = {dn11[0] - dn00[0], dn11[1] - dn00[1], dn11[2] - dn00[2]};
    double dd2[3] = {dn10[0] - dn01[0], dn10[1] - dn01[1], dn10[2] - dn01[2]};

    // d_raw = dd1 x d2 + d1 x dd2
    double d_raw[3];
    d_raw[0] =
        (dd1[1] * d2[2] - dd1[2] * d2[1]) + (d1[1] * dd2[2] - d1[2] * dd2[1]);
    d_raw[1] =
        (dd1[2] * d2[0] - dd1[0] * d2[2]) + (d1[2] * dd2[0] - d1[0] * dd2[2]);
    d_raw[2] =
        (dd1[0] * d2[1] - dd1[1] * d2[0]) + (d1[0] * dd2[1] - d1[1] * dd2[0]);

    // d(normal) = (d_raw - normal * (normal . d_raw)) / mag
    double n_dot_dr =
        normal[0] * d_raw[0] + normal[1] * d_raw[1] + normal[2] * d_raw[2];
    d_normal[0] = (d_raw[0] - normal[0] * n_dot_dr) * inv_mag;
    d_normal[1] = (d_raw[1] - normal[1] * n_dot_dr) * inv_mag;
    d_normal[2] = (d_raw[2] - normal[2] * n_dot_dr) * inv_mag;
}

/// Compute core radius for a vortex, matching solver and analytical_aic logic.
inline double compute_core_radius(double strip_width, bool same_component,
                                  double chord, double core_chord_fraction,
                                  double core_width_fraction) {
    if (same_component) {
        return 0.0001 * strip_width;
    }
    double chord_core = core_chord_fraction * chord;
    double width_core = core_width_fraction * strip_width;
    return (chord_core > width_core) ? chord_core : width_core;
}

inline double compute_strip_width(const double* rv1, const double* rv2) {
    double dy = rv2[1] - rv1[1];
    double dz = rv2[2] - rv1[2];
    return std::sqrt(dy * dy + dz * dz);
}

/// Compute JVP of horseshoe velocity w.r.t. one input parameter using 3 VJPs.
///
/// Given a VJP function vjp(seed, result) that computes result = J^T * seed,
/// compute the JVP: delta_V = J * delta_p.
///
/// Uses: (J * dp)_m = (J^T * e_m) . dp = vjp(e_m) . dp
template <typename VJPFunc>
void horseshoe_jvp_via_3vjps(VJPFunc vjp_func, const double* delta_p,
                             double* delta_v) {
    double seed[3], vjp_result[3];
    for (int m = 0; m < 3; ++m) {
        seed[0] = seed[1] = seed[2] = 0.0;
        seed[m] = 1.0;
        vjp_func(seed, vjp_result);
        delta_v[m] = vjp_result[0] * delta_p[0] + vjp_result[1] * delta_p[1] +
                     vjp_result[2] * delta_p[2];
    }
}

}  // namespace

// ============================================================================
// PropagateNodePerturbationsToVortex
// ============================================================================

VortexGeometryPerturbation PropagateNodePerturbationsToVortex(
    const SurfaceMesh& surface_mesh, Index first_vortex, Index nx, Index ny,
    const MeshNodeGradients& delta_nodes) {
    Index n_total = nx * ny;
    VortexGeometryPerturbation result;
    result.Resize(first_vortex + n_total);

    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            Index v = first_vortex + j * nx + i;

            // Extract node perturbations for this panel's 4 corners
            double dn00[3], dn10[3], dn01[3], dn11[3];
            for (int k = 0; k < 3; ++k) {
                dn00[k] = delta_nodes(i, j, k);
                dn10[k] = delta_nodes(i + 1, j, k);
                dn01[k] = delta_nodes(i, j + 1, k);
                dn11[k] = delta_nodes(i + 1, j + 1, k);
            }

            // rv1 = 0.75*N(i,j) + 0.25*N(i+1,j)
            for (int k = 0; k < 3; ++k) {
                result.d_rv1[3 * v + k] = 0.75 * dn00[k] + 0.25 * dn10[k];
            }

            // rv2 = 0.75*N(i,j+1) + 0.25*N(i+1,j+1)
            for (int k = 0; k < 3; ++k) {
                result.d_rv2[3 * v + k] = 0.75 * dn01[k] + 0.25 * dn11[k];
            }

            // rc = 0.5*(0.25*N(i,j) + 0.75*N(i+1,j))
            //    + 0.5*(0.25*N(i,j+1) + 0.75*N(i+1,j+1))
            for (int k = 0; k < 3; ++k) {
                result.d_rc[3 * v + k] = 0.125 * dn00[k] + 0.375 * dn10[k] +
                                         0.125 * dn01[k] + 0.375 * dn11[k];
            }

            // Normal: forward through normalize(cross(n11-n00, n10-n01))
            Vec3 p00 = surface_mesh.GetNode(i, j);
            Vec3 p10 = surface_mesh.GetNode(i + 1, j);
            Vec3 p01 = surface_mesh.GetNode(i, j + 1);
            Vec3 p11 = surface_mesh.GetNode(i + 1, j + 1);

            double n00c[3] = {p00.x(), p00.y(), p00.z()};
            double n10c[3] = {p10.x(), p10.y(), p10.z()};
            double n01c[3] = {p01.x(), p01.y(), p01.z()};
            double n11c[3] = {p11.x(), p11.y(), p11.z()};

            ForwardNormalPerturbation(n00c, n10c, n01c, n11c, dn00, dn10, dn01,
                                      dn11, &result.d_normal[3 * v]);
        }
    }

    return result;
}

// ============================================================================
// ComputeResidualGeometryJVP
// ============================================================================

VecX ComputeResidualGeometryJVP(const VortexLatticeSolver& solver,
                                const SolutionState& solution,
                                const FlightCondition& condition,
                                const VortexGeometryPerturbation& delta_geom) {
    const VortexLattice& mesh = solver.GetMesh();
    const Index n = mesh.GetNumVortices();

    VecX delta_R = VecX::Zero(n);

    // Extract mesh arrays
    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta_pg =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    // Freestream direction
    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    // Part 1: Row gradient for rc and normal contributions.
    // For each row i, the row gradient with weights=Gamma gives:
    //   d_rc = d(AIC_row . Gamma)/d(rc_i)
    //   d_normal_aic = d(AIC_row . Gamma)/d(normal_i) = induced velocity at i
    for (Index i = 0; i < n; ++i) {
        double d_rc_out[3] = {0, 0, 0};
        double d_normal_aic[3] = {0, 0, 0};

        analytical_aic_row_gradient_rc(
            arrays.rv1.data(), arrays.rv2.data(), &arrays.rc[3 * i],
            &arrays.normal[3 * i], arrays.chord.data(),
            arrays.component_idx.data(), static_cast<int>(n),
            static_cast<int>(i), arrays.component_idx[i], core_chord_frac,
            core_width_frac, beta_pg, solution.circulation.data(), d_rc_out,
            d_normal_aic);

        // Forward mode: delta(AIC*Gamma)[i] from rc perturbation
        delta_R(i) += d_rc_out[0] * delta_geom.d_rc[3 * i + 0] +
                      d_rc_out[1] * delta_geom.d_rc[3 * i + 1] +
                      d_rc_out[2] * delta_geom.d_rc[3 * i + 2];

        // Normal contribution: delta(AIC*Gamma)[i] from normal perturbation
        // minus delta(b)[i] from normal perturbation (b = -V_inf . normal)
        // Combined: (d_normal_aic + V_inf) . delta_normal
        // (Note: b = -V_inf . n, so delta_b = -V_inf . delta_n,
        //  and delta_R = delta(AIC*Gamma) - delta_b = d_normal_aic . dn + V_inf
        //  . dn)
        delta_R(i) +=
            (d_normal_aic[0] + v_inf[0]) * delta_geom.d_normal[3 * i + 0] +
            (d_normal_aic[1] + v_inf[1]) * delta_geom.d_normal[3 * i + 1] +
            (d_normal_aic[2] + v_inf[2]) * delta_geom.d_normal[3 * i + 2];
    }

    // Part 2: Column contributions from rv1/rv2 perturbations.
    // For each (i,j) pair, the perturbation of AIC[i,j] from rv1_j, rv2_j
    // changes uses the horseshoe velocity VJP seeded with normal_i.
    for (Index j = 0; j < n; ++j) {
        double gamma_j = solution.circulation(j);
        if (std::abs(gamma_j) < 1e-30)
            continue;

        const double* rv1_j = &arrays.rv1[3 * j];
        const double* rv2_j = &arrays.rv2[3 * j];
        double strip_width_j = compute_strip_width(rv1_j, rv2_j);
        int component_j = arrays.component_idx[j];

        // Skip if no perturbation for this vortex
        double drv1_norm = std::abs(delta_geom.d_rv1[3 * j]) +
                           std::abs(delta_geom.d_rv1[3 * j + 1]) +
                           std::abs(delta_geom.d_rv1[3 * j + 2]);
        double drv2_norm = std::abs(delta_geom.d_rv2[3 * j]) +
                           std::abs(delta_geom.d_rv2[3 * j + 1]) +
                           std::abs(delta_geom.d_rv2[3 * j + 2]);
        if (drv1_norm + drv2_norm < 1e-30)
            continue;

        for (Index i = 0; i < n; ++i) {
            const double* rc_i = &arrays.rc[3 * i];
            const double* normal_i = &arrays.normal[3 * i];

            bool same_component = (component_j == arrays.component_idx[i]);
            double core_radius = compute_core_radius(
                strip_width_j, same_component, arrays.chord[j], core_chord_frac,
                core_width_frac);
            bool include_bound = (i != j);

            // VJP of horseshoe velocity seeded with normal gives gradient
            // of AIC element w.r.t. rv1, rv2
            double vjp_rv1[3], vjp_rv2[3];
            horseshoe_velocity_vjp_bound_left(rc_i, rv1_j, rv2_j, beta_pg,
                                              core_radius, include_bound,
                                              normal_i, vjp_rv1);
            horseshoe_velocity_vjp_bound_right(rc_i, rv1_j, rv2_j, beta_pg,
                                               core_radius, include_bound,
                                               normal_i, vjp_rv2);

            // Forward: delta(AIC[i,j]) from rv1_j, rv2_j perturbation
            double dA_ij = vjp_rv1[0] * delta_geom.d_rv1[3 * j + 0] +
                           vjp_rv1[1] * delta_geom.d_rv1[3 * j + 1] +
                           vjp_rv1[2] * delta_geom.d_rv1[3 * j + 2] +
                           vjp_rv2[0] * delta_geom.d_rv2[3 * j + 0] +
                           vjp_rv2[1] * delta_geom.d_rv2[3 * j + 1] +
                           vjp_rv2[2] * delta_geom.d_rv2[3 * j + 2];

            delta_R(i) += gamma_j * dA_ij;
        }
    }

    return delta_R;
}

// ============================================================================
// ComputePerVortexForceTotalJVP
// ============================================================================

VecX ComputePerVortexForceTotalJVP(const VortexLatticeSolver& solver,
                                   const SolutionState& solution,
                                   const FlightCondition& condition,
                                   const VortexGeometryPerturbation& delta_geom,
                                   const VecX& delta_gamma) {
    const VortexLattice& mesh = solver.GetMesh();
    const MatX& vim = solver.GetVelocityInfluenceMatrix();
    const Index n = mesh.GetNumVortices();

    VecX delta_F = VecX::Zero(3 * n);

    // Extract mesh arrays
    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta_pg =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    for (Index i = 0; i < n; ++i) {
        double gamma_i = solution.circulation(i);

        // Current effective velocity and bound segment
        double v_eff[3] = {v_inf[0] + solution.induced_velocity(i, 0),
                           v_inf[1] + solution.induced_velocity(i, 1),
                           v_inf[2] + solution.induced_velocity(i, 2)};

        double dl[3] = {arrays.rv2[3 * i + 0] - arrays.rv1[3 * i + 0],
                        arrays.rv2[3 * i + 1] - arrays.rv1[3 * i + 1],
                        arrays.rv2[3 * i + 2] - arrays.rv1[3 * i + 2]};

        // ============================================================
        // Term 1: 2*delta_Gamma_i * (V_eff_i x dl_i)
        // ============================================================
        double v_cross_dl[3] = {v_eff[1] * dl[2] - v_eff[2] * dl[1],
                                v_eff[2] * dl[0] - v_eff[0] * dl[2],
                                v_eff[0] * dl[1] - v_eff[1] * dl[0]};

        for (int k = 0; k < 3; ++k) {
            delta_F[3 * i + k] += 2.0 * delta_gamma(i) * v_cross_dl[k];
        }

        // ============================================================
        // Term 2: 2*Gamma_i * (V_eff_i x delta_dl_i)
        // where delta_dl = delta_rv2 - delta_rv1
        // ============================================================
        double d_dl[3] = {
            delta_geom.d_rv2[3 * i + 0] - delta_geom.d_rv1[3 * i + 0],
            delta_geom.d_rv2[3 * i + 1] - delta_geom.d_rv1[3 * i + 1],
            delta_geom.d_rv2[3 * i + 2] - delta_geom.d_rv1[3 * i + 2]};

        double v_cross_ddl[3] = {v_eff[1] * d_dl[2] - v_eff[2] * d_dl[1],
                                 v_eff[2] * d_dl[0] - v_eff[0] * d_dl[2],
                                 v_eff[0] * d_dl[1] - v_eff[1] * d_dl[0]};

        for (int k = 0; k < 3; ++k) {
            delta_F[3 * i + k] += 2.0 * gamma_i * v_cross_ddl[k];
        }

        // ============================================================
        // Term 3: 2*Gamma_i * (delta_V_ind_i x dl_i)
        // delta_V_ind = delta_VIM*Gamma + VIM*delta_Gamma
        // ============================================================

        // Part 3a: VIM * delta_Gamma contribution
        double dv_from_dgamma[3] = {0.0, 0.0, 0.0};
        for (Index jj = 0; jj < n; ++jj) {
            for (int k = 0; k < 3; ++k) {
                dv_from_dgamma[k] += vim(3 * i + k, jj) * delta_gamma(jj);
            }
        }

        // Part 3b: delta_VIM * Gamma contribution (VIM geometry sensitivity)
        // For each column j, the horseshoe velocity from vortex j at rc_i
        // changes when rc_i, rv1_j, rv2_j are perturbed.
        double dv_from_geom[3] = {0.0, 0.0, 0.0};

        for (Index j = 0; j < n; ++j) {
            double gamma_j = solution.circulation(j);
            if (std::abs(gamma_j) < 1e-30)
                continue;

            const double* rc_i_ptr = &arrays.rc[3 * i];
            const double* rv1_j_ptr = &arrays.rv1[3 * j];
            const double* rv2_j_ptr = &arrays.rv2[3 * j];

            bool same_component =
                (arrays.component_idx[j] == arrays.component_idx[i]);
            double sw = compute_strip_width(rv1_j_ptr, rv2_j_ptr);
            double core_radius =
                compute_core_radius(sw, same_component, arrays.chord[j],
                                    core_chord_frac, core_width_frac);
            bool include_bound = (i != j);

            // Compute JVP of horseshoe velocity via 3 VJPs per derivative type.
            // For field_point (rc) perturbation:
            double dv_rc[3] = {0, 0, 0};
            if (std::abs(delta_geom.d_rc[3 * i]) +
                    std::abs(delta_geom.d_rc[3 * i + 1]) +
                    std::abs(delta_geom.d_rc[3 * i + 2]) >
                1e-30) {
                horseshoe_jvp_via_3vjps(
                    [&](const double* seed, double* res) {
                        horseshoe_velocity_vjp_field_point(
                            rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg,
                            core_radius, include_bound, seed, res);
                    },
                    &delta_geom.d_rc[3 * i], dv_rc);
            }

            // For bound_left (rv1_j) perturbation:
            double dv_rv1[3] = {0, 0, 0};
            if (std::abs(delta_geom.d_rv1[3 * j]) +
                    std::abs(delta_geom.d_rv1[3 * j + 1]) +
                    std::abs(delta_geom.d_rv1[3 * j + 2]) >
                1e-30) {
                horseshoe_jvp_via_3vjps(
                    [&](const double* seed, double* res) {
                        horseshoe_velocity_vjp_bound_left(
                            rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg,
                            core_radius, include_bound, seed, res);
                    },
                    &delta_geom.d_rv1[3 * j], dv_rv1);
            }

            // For bound_right (rv2_j) perturbation:
            double dv_rv2[3] = {0, 0, 0};
            if (std::abs(delta_geom.d_rv2[3 * j]) +
                    std::abs(delta_geom.d_rv2[3 * j + 1]) +
                    std::abs(delta_geom.d_rv2[3 * j + 2]) >
                1e-30) {
                horseshoe_jvp_via_3vjps(
                    [&](const double* seed, double* res) {
                        horseshoe_velocity_vjp_bound_right(
                            rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg,
                            core_radius, include_bound, seed, res);
                    },
                    &delta_geom.d_rv2[3 * j], dv_rv2);
            }

            // Accumulate: delta_VIM[i,j] * Gamma_j
            for (int k = 0; k < 3; ++k) {
                dv_from_geom[k] += gamma_j * (dv_rc[k] + dv_rv1[k] + dv_rv2[k]);
            }
        }

        // Total induced velocity perturbation
        double dv_ind[3];
        for (int k = 0; k < 3; ++k) {
            dv_ind[k] = dv_from_dgamma[k] + dv_from_geom[k];
        }

        // delta_V_ind x dl
        double dv_cross_dl[3] = {dv_ind[1] * dl[2] - dv_ind[2] * dl[1],
                                 dv_ind[2] * dl[0] - dv_ind[0] * dl[2],
                                 dv_ind[0] * dl[1] - dv_ind[1] * dl[0]};

        for (int k = 0; k < 3; ++k) {
            delta_F[3 * i + k] += 2.0 * gamma_i * dv_cross_dl[k];
        }
    }

    return delta_F;
}

// ============================================================================
// ComputeResidualGeometryJVPBatch
// ============================================================================

MatX ComputeResidualGeometryJVPBatch(
    const VortexLatticeSolver& solver, const SolutionState& solution,
    const FlightCondition& condition,
    const VortexGeometryPerturbationBatch& batch) {
    const VortexLattice& mesh = solver.GetMesh();
    const Index n = mesh.GetNumVortices();
    const Index kk = batch.NumDirections();

    MatX delta_R = MatX::Zero(n, kk);

    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta_pg =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    // Row-parallel sweep: each row i evaluates its kernels once and
    // contracts them with every perturbation direction
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
    for (Index i = 0; i < n; ++i) {
        // Part 1: rc and normal contributions (row gradient, weights=Gamma)
        double d_rc_out[3] = {0, 0, 0};
        double d_normal_aic[3] = {0, 0, 0};

        analytical_aic_row_gradient_rc(
            arrays.rv1.data(), arrays.rv2.data(), &arrays.rc[3 * i],
            &arrays.normal[3 * i], arrays.chord.data(),
            arrays.component_idx.data(), static_cast<int>(n),
            static_cast<int>(i), arrays.component_idx[i], core_chord_frac,
            core_width_frac, beta_pg, solution.circulation.data(), d_rc_out,
            d_normal_aic);

        Eigen::RowVector3d w_rc(d_rc_out[0], d_rc_out[1], d_rc_out[2]);
        Eigen::RowVector3d w_n(d_normal_aic[0] + v_inf[0],
                               d_normal_aic[1] + v_inf[1],
                               d_normal_aic[2] + v_inf[2]);

        delta_R.row(i) += w_rc * batch.d_rc.middleRows(3 * i, 3);
        delta_R.row(i) += w_n * batch.d_normal.middleRows(3 * i, 3);

        // Part 2: rv1/rv2 column contributions. The VJPs seeded with
        // normal_i give d(AIC[i,j])/d(rv1_j, rv2_j); evaluate once per
        // pair and contract with all directions.
        const double* rc_i = &arrays.rc[3 * i];
        const double* normal_i = &arrays.normal[3 * i];

        for (Index j = 0; j < n; ++j) {
            double gamma_j = solution.circulation(j);
            if (std::abs(gamma_j) < 1e-30)
                continue;

            const double* rv1_j = &arrays.rv1[3 * j];
            const double* rv2_j = &arrays.rv2[3 * j];
            double strip_width_j = compute_strip_width(rv1_j, rv2_j);
            bool same_component =
                (arrays.component_idx[j] == arrays.component_idx[i]);
            double core_radius = compute_core_radius(
                strip_width_j, same_component, arrays.chord[j], core_chord_frac,
                core_width_frac);
            bool include_bound = (i != j);

            double vjp_rv1[3], vjp_rv2[3];
            horseshoe_velocity_vjp_bound_left(rc_i, rv1_j, rv2_j, beta_pg,
                                              core_radius, include_bound,
                                              normal_i, vjp_rv1);
            horseshoe_velocity_vjp_bound_right(rc_i, rv1_j, rv2_j, beta_pg,
                                               core_radius, include_bound,
                                               normal_i, vjp_rv2);

            Eigen::RowVector3d w1(gamma_j * vjp_rv1[0], gamma_j * vjp_rv1[1],
                                  gamma_j * vjp_rv1[2]);
            Eigen::RowVector3d w2(gamma_j * vjp_rv2[0], gamma_j * vjp_rv2[1],
                                  gamma_j * vjp_rv2[2]);

            delta_R.row(i) += w1 * batch.d_rv1.middleRows(3 * j, 3);
            delta_R.row(i) += w2 * batch.d_rv2.middleRows(3 * j, 3);
        }
    }

    return delta_R;
}

// ============================================================================
// ComputePerVortexForceTotalJVPBatch
// ============================================================================

MatX ComputePerVortexForceTotalJVPBatch(
    const VortexLatticeSolver& solver, const SolutionState& solution,
    const FlightCondition& condition,
    const VortexGeometryPerturbationBatch& batch, const MatX& delta_gammas) {
    const VortexLattice& mesh = solver.GetMesh();
    const MatX& vim = solver.GetVelocityInfluenceMatrix();
    const Index n = mesh.GetNumVortices();
    const Index kk = batch.NumDirections();

    MatX delta_F = MatX::Zero(3 * n, kk);

    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta_pg =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    // Term 3a for all vortices and directions in one product:
    // delta_V_ind (from delta_Gamma) = VIM * delta_Gamma
    MatX dv_from_dgamma = vim * delta_gammas;  // (3N x K)

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8)
#endif
    for (Index i = 0; i < n; ++i) {
        double gamma_i = solution.circulation(i);

        Eigen::Vector3d v_eff(v_inf[0] + solution.induced_velocity(i, 0),
                              v_inf[1] + solution.induced_velocity(i, 1),
                              v_inf[2] + solution.induced_velocity(i, 2));
        Eigen::Vector3d dl(arrays.rv2[3 * i + 0] - arrays.rv1[3 * i + 0],
                           arrays.rv2[3 * i + 1] - arrays.rv1[3 * i + 1],
                           arrays.rv2[3 * i + 2] - arrays.rv1[3 * i + 2]);

        Eigen::Matrix3d skew_veff;
        skew_veff << 0, -v_eff.z(), v_eff.y(), v_eff.z(), 0, -v_eff.x(),
            -v_eff.y(), v_eff.x(), 0;
        Eigen::Matrix3d skew_dl;
        skew_dl << 0, -dl.z(), dl.y(), dl.z(), 0, -dl.x(), -dl.y(), dl.x(), 0;

        // Term 1: 2 * delta_Gamma_i * (V_eff x dl) for every direction
        Eigen::Vector3d v_cross_dl = v_eff.cross(dl);
        delta_F.middleRows(3 * i, 3).noalias() +=
            2.0 * v_cross_dl * delta_gammas.row(i);

        // Term 2: 2 * Gamma_i * (V_eff x delta_dl)
        delta_F.middleRows(3 * i, 3).noalias() +=
            2.0 * gamma_i * skew_veff *
            (batch.d_rv2.middleRows(3 * i, 3) -
             batch.d_rv1.middleRows(3 * i, 3));

        // Term 3: 2 * Gamma_i * (delta_V_ind x dl), with
        // delta_V_ind = VIM*delta_Gamma (precomputed) + delta_VIM*Gamma.
        // The geometry part evaluates the pair kernels once, assembling
        // 3x3 Jacobians applied to all directions. Contributions through
        // rc_i are shared across j and accumulate into one 3x3.
        Eigen::MatrixXd dv_geom = Eigen::MatrixXd::Zero(3, kk);
        Eigen::Matrix3d J_rc_sum = Eigen::Matrix3d::Zero();

        const double* rc_i_ptr = &arrays.rc[3 * i];

        for (Index j = 0; j < n; ++j) {
            double gamma_j = solution.circulation(j);
            if (std::abs(gamma_j) < 1e-30)
                continue;

            const double* rv1_j_ptr = &arrays.rv1[3 * j];
            const double* rv2_j_ptr = &arrays.rv2[3 * j];

            bool same_component =
                (arrays.component_idx[j] == arrays.component_idx[i]);
            double sw = compute_strip_width(rv1_j_ptr, rv2_j_ptr);
            double core_radius =
                compute_core_radius(sw, same_component, arrays.chord[j],
                                    core_chord_frac, core_width_frac);
            bool include_bound = (i != j);

            // 3x3 Jacobians of the horseshoe velocity: row m = VJP(e_m)
            double seed[3], res[3];
            Eigen::Matrix3d J_rc, J_rv1, J_rv2;
            for (int m = 0; m < 3; ++m) {
                seed[0] = seed[1] = seed[2] = 0.0;
                seed[m] = 1.0;
                horseshoe_velocity_vjp_field_point(
                    rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg, core_radius,
                    include_bound, seed, res);
                J_rc.row(m) << res[0], res[1], res[2];
                horseshoe_velocity_vjp_bound_left(
                    rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg, core_radius,
                    include_bound, seed, res);
                J_rv1.row(m) << res[0], res[1], res[2];
                horseshoe_velocity_vjp_bound_right(
                    rc_i_ptr, rv1_j_ptr, rv2_j_ptr, beta_pg, core_radius,
                    include_bound, seed, res);
                J_rv2.row(m) << res[0], res[1], res[2];
            }

            J_rc_sum.noalias() += gamma_j * J_rc;
            dv_geom.noalias() +=
                gamma_j * J_rv1 * batch.d_rv1.middleRows(3 * j, 3);
            dv_geom.noalias() +=
                gamma_j * J_rv2 * batch.d_rv2.middleRows(3 * j, 3);
        }

        dv_geom.noalias() += J_rc_sum * batch.d_rc.middleRows(3 * i, 3);

        // (delta_V_ind x dl) = -dl x delta_V_ind = -skew(dl) * delta_V_ind
        delta_F.middleRows(3 * i, 3).noalias() +=
            2.0 * gamma_i * (-skew_dl) *
            (dv_from_dgamma.middleRows(3 * i, 3) + dv_geom);
    }

    return delta_F;
}

}  // namespace pegasus::gradients
