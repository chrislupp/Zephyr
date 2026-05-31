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

#include "pegasus/viscous/coupler.h"

#include "pegasus/solver/vortex_lattice.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

namespace pegasus::viscous {

namespace {

/// Compute α_eff (radians) at a strip from the VLM solution.
///
/// Defines effective angle as the angle between the strip's chord line and
/// the local total velocity (freestream + induced). For a planar wing this
/// is just  atan2(v_total.z, v_total.x) − strip.twist_angle. For swept or
/// dihedral wings additional projection onto the strip's chord plane would
/// be required (a v1 limitation noted in the coupler header).
Scalar ComputeStripEffectiveAlpha(const VortexLattice& mesh,
                                  const SolutionState& solution,
                                  const FlightCondition& condition,
                                  Index strip_idx) {
    const Strip& strip = mesh.GetStrip(strip_idx);
    Index iv = strip.first_vortex_index;

    Vec3 v_inf = VortexLatticeSolver::ComputeFreestreamDirection(condition);
    Vec3 v_ind(solution.induced_velocity(iv, 0),
               solution.induced_velocity(iv, 1),
               solution.induced_velocity(iv, 2));
    Vec3 v_total = v_inf + v_ind;

    Scalar alpha_body = std::atan2(v_total.z(), v_total.x());
    return alpha_body - strip.twist_angle;
}

}  // namespace

ViscousVlmCoupler::ViscousVlmCoupler(VortexLatticeSolver& solver,
                                     const FlightCondition& flight,
                                     const CouplerOptions& options)
    : solver_(solver), flight_(flight), options_(options) {
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_strips = mesh.GetNumStrips();

    baseline_twist_.resize(n_strips);
    delta_alpha_ = VecX::Zero(n_strips);
    cd_2d_ = VecX::Zero(n_strips);

    strip_polars_.reserve(static_cast<std::size_t>(n_strips));
    for (Index i = 0; i < n_strips; ++i) {
        baseline_twist_(i) = mesh.GetStrip(i).twist_angle;
        strip_polars_.emplace_back(options_.airfoil_naca,
                                   options_.scirocco_panels);
    }
}

void ViscousVlmCoupler::ApplyTwistsAndRefresh() {
    VortexLattice& mesh = solver_.GetMeshMutable();
    const Index n = mesh.GetNumStrips();
    for (Index i = 0; i < n; ++i) {
        mesh.GetStrip(i).twist_angle = baseline_twist_(i) + delta_alpha_(i);
    }
    mesh.RecalculateNormals();
    solver_.RebuildAIC();
}

void ViscousVlmCoupler::EvaluateVlm(VecX& cl_vlm, VecX& alpha_eff) {
    ApplyTwistsAndRefresh();
    SolutionState solution = solver_.Solve(flight_);
    last_results_ = force_calc_.ComputeForces(
        solver_.GetReferenceQuantities(), solver_.GetMesh(), solution, flight_);

    const Index n = static_cast<Index>(strip_polars_.size());
    cl_vlm.resize(n);
    alpha_eff.resize(n);

    // Flatten the surface→strip nested arrangement into a flat strip index
    // matching mesh.GetStrip(i). Build by walking surface_loads in mesh order.
    const VortexLattice& mesh = solver_.GetMesh();
    for (Index strip_idx = 0; strip_idx < n; ++strip_idx) {
        alpha_eff(strip_idx) =
            ComputeStripEffectiveAlpha(mesh, solution, flight_, strip_idx);

        // Find this strip's StripLoads inside last_results_.surface_loads.
        const Strip& strip = mesh.GetStrip(strip_idx);
        Index surf = strip.surface_index;
        const auto& surf_loads = last_results_.surface_loads[surf];
        // Strips inside a surface are appended in mesh order; find by
        // strip_index.
        Scalar cl_local = 0.0;
        for (const auto& s : surf_loads.strip_loads) {
            if (s.strip_index == strip_idx) {
                cl_local = s.cl;
                break;
            }
        }
        cl_vlm(strip_idx) = cl_local;
    }
}

VecX ViscousVlmCoupler::ComputeResidual(const VecX& alpha_eff,
                                        const VecX& cl_vlm) {
    const Index n = alpha_eff.size();
    VecX r(n);
    cd_2d_.resize(n);
    for (Index i = 0; i < n; ++i) {
        bool ok = strip_polars_[static_cast<std::size_t>(i)].SolveFixedAlpha(
            alpha_eff(i), options_.reynolds, options_.mach);
        if (!ok) {
            // Fallback: use thin-airfoil estimate if Scirocco fails. This
            // keeps the Newton loop alive so the user gets some answer back.
            Scalar cl2 = 2.0 * kPi * alpha_eff(i);
            r(i) = cl_vlm(i) - cl2;
            cd_2d_(i) = 0.0;
        } else {
            Scalar cl2 = strip_polars_[static_cast<std::size_t>(i)].GetCl();
            cd_2d_(i) = strip_polars_[static_cast<std::size_t>(i)].GetCd();
            r(i) = cl_vlm(i) - cl2;
        }
    }
    return r;
}

MatX ViscousVlmCoupler::BuildJacobianFD(const VecX& base_residual) {
    const Index n = base_residual.size();
    MatX J(n, n);

    VecX cl_pert(n), alpha_pert(n);
    const VecX delta_alpha_save = delta_alpha_;

    // For each column j, perturb Δα_j by h, re-solve VLM (Scirocco terms
    // come along), record new residual.
    for (Index j = 0; j < n; ++j) {
        delta_alpha_ = delta_alpha_save;
        delta_alpha_(j) += options_.fd_step_rad;

        EvaluateVlm(cl_pert, alpha_pert);
        VecX r_pert = ComputeResidual(alpha_pert, cl_pert);

        J.col(j) = (r_pert - base_residual) / options_.fd_step_rad;
    }
    // Restore baseline state.
    delta_alpha_ = delta_alpha_save;
    return J;
}

CouplerResult ViscousVlmCoupler::Solve() {
    CouplerResult out;
    const Index n = static_cast<Index>(strip_polars_.size());
    if (n == 0) {
        out.converged = true;
        return out;
    }

    VecX cl_vlm(n), alpha_eff(n);

    for (int iter = 0; iter < options_.max_iterations; ++iter) {
        EvaluateVlm(cl_vlm, alpha_eff);
        VecX r = ComputeResidual(alpha_eff, cl_vlm);
        Scalar res_norm = r.lpNorm<Eigen::Infinity>();
        out.iterations = iter + 1;
        out.final_residual = res_norm;

        if (res_norm < options_.tolerance) {
            out.converged = true;
            break;
        }

        // Build FD Jacobian (this re-evaluates VLM+Scirocco N times — the
        // most expensive part of v1).
        MatX J = BuildJacobianFD(r);

        // Newton step: solve J Δu = −r
        Eigen::FullPivLU<MatX> lu(J);
        VecX du = lu.solve(-r);

        // Apply with damping.
        delta_alpha_ += options_.damping * du;
    }

    // Final evaluation at the converged Δα to refresh results.
    EvaluateVlm(cl_vlm, alpha_eff);
    (void)ComputeResidual(alpha_eff, cl_vlm);

    // Add viscous drag to the aerodynamic results. CD_viscous = sum over
    // strips of (Cd_2D,i * strip_area_i) / reference_area.
    const VortexLattice& mesh = solver_.GetMesh();
    Scalar cdv = 0.0;
    Scalar ref_area = last_results_.reference_area;
    for (Index i = 0; i < n; ++i) {
        const Strip& strip = mesh.GetStrip(i);
        Scalar strip_area = strip.chord_center * strip.width;
        cdv += cd_2d_(i) * strip_area;
    }
    cdv /= std::max(ref_area, Scalar(1e-30));
    last_results_.CDv = cdv;
    // Total CD with viscous correction.
    last_results_.CD = last_results_.CDi_nearfield + cdv;

    out.results = last_results_;
    out.delta_alpha = delta_alpha_;
    return out;
}

CoupledGradients ViscousVlmCoupler::ComputeAlphaDerivatives() {
    CoupledGradients g;
    const Index n = static_cast<Index>(strip_polars_.size());

    // Base evaluation at converged (u*, α∞)
    VecX cl_base(n), alpha_eff_base(n);
    EvaluateVlm(cl_base, alpha_eff_base);
    VecX r_base = ComputeResidual(alpha_eff_base, cl_base);
    AerodynamicResults base_results = last_results_;
    // Capture base CL/CD with viscous correction.
    auto AccumulateCDv = [&](const VecX& cd2d) {
        Scalar cdv = 0.0;
        const VortexLattice& mesh = solver_.GetMesh();
        for (Index i = 0; i < n; ++i) {
            const Strip& strip = mesh.GetStrip(i);
            cdv += cd2d(i) * strip.chord_center * strip.width;
        }
        return cdv / std::max(base_results.reference_area, Scalar(1e-30));
    };
    const VecX cd2d_base = cd_2d_;
    const Scalar CL_base = base_results.CL;
    const Scalar CD_base =
        base_results.CDi_nearfield + AccumulateCDv(cd2d_base);
    g.CL = CL_base;
    g.CD = CD_base;

    // -----------------------------------------------------------------
    // ∂F/∂u and J = ∂R/∂u via finite differences on Δα (shared sweep).
    // For each strip j perturb Δα_j by h, re-evaluate everything, capture
    // both r_pert (residual column) and F_pert (CL, CD with viscous
    // correction). Net cost: N evaluations.
    // -----------------------------------------------------------------
    const Scalar h = options_.fd_step_rad;
    MatX J(n, n);  // ∂R/∂u
    VecX dCL_du = VecX::Zero(n);
    VecX dCD_du = VecX::Zero(n);

    const VecX delta_alpha_save = delta_alpha_;
    VecX cl_pert(n), alpha_pert(n);

    for (Index j = 0; j < n; ++j) {
        delta_alpha_ = delta_alpha_save;
        delta_alpha_(j) += h;

        EvaluateVlm(cl_pert, alpha_pert);
        VecX r_pert = ComputeResidual(alpha_pert, cl_pert);
        J.col(j) = (r_pert - r_base) / h;

        const Scalar CL_pert = last_results_.CL;
        const Scalar CD_pert =
            last_results_.CDi_nearfield + AccumulateCDv(cd_2d_);
        dCL_du(j) = (CL_pert - CL_base) / h;
        dCD_du(j) = (CD_pert - CD_base) / h;
    }
    delta_alpha_ = delta_alpha_save;

    // -----------------------------------------------------------------
    // ∂F/∂α∞ and ∂R/∂α∞ via FD on α∞, holding u fixed.
    // -----------------------------------------------------------------
    const Scalar h_alpha = options_.fd_step_rad;
    const FlightCondition flight_save = flight_;
    flight_.alpha += h_alpha;
    EvaluateVlm(cl_pert, alpha_pert);
    VecX r_alpha = ComputeResidual(alpha_pert, cl_pert);
    const Scalar CL_alpha_pert = last_results_.CL;
    const Scalar CD_alpha_pert =
        last_results_.CDi_nearfield + AccumulateCDv(cd_2d_);
    flight_ = flight_save;

    VecX dR_dalpha = (r_alpha - r_base) / h_alpha;
    Scalar dCL_dalpha_partial = (CL_alpha_pert - CL_base) / h_alpha;
    Scalar dCD_dalpha_partial = (CD_alpha_pert - CD_base) / h_alpha;

    // -----------------------------------------------------------------
    // Forward sensitivity:  J · (du/dα) = −∂R/∂α
    // -----------------------------------------------------------------
    Eigen::FullPivLU<MatX> lu(J);
    VecX du_dalpha = lu.solve(-dR_dalpha);

    // -----------------------------------------------------------------
    // Total derivative:  dF/dα = ∂F/∂α + ∂F/∂u · (du/dα)
    // -----------------------------------------------------------------
    g.dCL_dalpha = dCL_dalpha_partial + dCL_du.dot(du_dalpha);
    g.dCD_dalpha = dCD_dalpha_partial + dCD_du.dot(du_dalpha);

    // Restore the converged base state so subsequent calls see a clean
    // coupler.
    EvaluateVlm(cl_base, alpha_eff_base);
    (void)ComputeResidual(alpha_eff_base, cl_base);
    last_results_ = base_results;
    cd_2d_ = cd2d_base;

    return g;
}

}  // namespace pegasus::viscous
