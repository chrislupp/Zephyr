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

#include "pegasus/gradients/analytical_forces.h"
#include "pegasus/solver/vortex_lattice.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

namespace pegasus::viscous {

namespace {

/// Compute α_eff (radians) at a strip from a VLM solution.
///
/// α_eff is the angle Scirocco is called at to match the per-strip Cl. In
/// Pegasus's per-strip Cl convention (lift is `cf.z·cos α∞ - cf.x·sin α∞`,
/// projected onto the *stability* axis defined by the global freestream),
/// the consistent definition is the body-frame angle of the total local
/// velocity, *minus* the strip's geometric twist. This makes Newton on the
/// residual `Cl_VLM − Cl_2D(α_eff)` well-conditioned (J is diagonal-
/// dominant) and converges to a physical state.
///
/// A purely geometric "+ twist" formula is more intuitive but gives a
/// poorly-conditioned Jacobian (3D lift slope cancels 2π in J_ii) and
/// Newton walks to a non-physical lift-free fixed point on small wings.
/// The "− twist" sign matches the implicit stability-axis projection
/// already baked into the per-strip Cl.
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

/// d/d(ang) of the normalized control-point normal, evaluated at the current
/// ang_cp = strip.twist_angle − atan(camber_slope). Since chord_dir is the
/// only twist-dependent factor in normal = normalize(chord_dir × bound_dir),
/// this is purely closed-form.
Vec3 NormalDerivativeWrtTwist(const HorseshoeVortex& vor, Scalar twist_angle,
                              Scalar trefftz_ny, Scalar trefftz_nz) {
    Scalar ang = twist_angle - std::atan(vor.camber_slope_at_control_point);
    Scalar ca = std::cos(ang);
    Scalar sa = std::sin(ang);

    Vec3 chord_dir(ca, -sa * trefftz_ny, -sa * trefftz_nz);
    Vec3 d_chord(-sa, -ca * trefftz_ny, -ca * trefftz_nz);

    Vec3 bound_raw = vor.bound_vortex_right - vor.bound_vortex_left;
    Scalar bound_norm = bound_raw.norm();
    if (bound_norm < Scalar(1e-30)) {
        return Vec3::Zero();
    }
    Vec3 bound_dir = bound_raw / bound_norm;

    Vec3 n_unnorm = chord_dir.cross(bound_dir);
    Scalar n_mag = n_unnorm.norm();
    if (n_mag < Scalar(1e-30)) {
        return Vec3::Zero();
    }
    Vec3 n = n_unnorm / n_mag;
    Vec3 d_n_unnorm = d_chord.cross(bound_dir);
    return (d_n_unnorm - n * (n.dot(d_n_unnorm))) / n_mag;
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
    dcl2d_dalpha_ = VecX::Zero(n_strips);
    dcd2d_dalpha_ = VecX::Zero(n_strips);

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
    last_solution_ = solver_.Solve(flight_);
    last_results_ =
        force_calc_.ComputeForces(solver_.GetReferenceQuantities(),
                                  solver_.GetMesh(), last_solution_, flight_);

    const Index n = static_cast<Index>(strip_polars_.size());
    cl_vlm.resize(n);
    alpha_eff.resize(n);

    const VortexLattice& mesh = solver_.GetMesh();
    for (Index strip_idx = 0; strip_idx < n; ++strip_idx) {
        alpha_eff(strip_idx) = ComputeStripEffectiveAlpha(mesh, last_solution_,
                                                          flight_, strip_idx);

        const Strip& strip = mesh.GetStrip(strip_idx);
        const auto& surf_loads =
            last_results_.surface_loads[strip.surface_index];
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
    dcl2d_dalpha_.resize(n);
    dcd2d_dalpha_.resize(n);

    for (Index i = 0; i < n; ++i) {
        auto& polar = strip_polars_[static_cast<std::size_t>(i)];
        bool ok = polar.SolveFixedAlpha(alpha_eff(i), options_.reynolds,
                                        options_.mach);
        if (!ok) {
            // Fallback: use thin-airfoil estimate so the Newton loop survives
            // an occasional Scirocco failure. dCl_2D/dα defaults to 2π.
            Scalar cl2 = 2.0 * kPi * alpha_eff(i);
            r(i) = cl_vlm(i) - cl2;
            cd_2d_(i) = 0.0;
            dcl2d_dalpha_(i) = 2.0 * kPi;
            dcd2d_dalpha_(i) = 0.0;
        } else {
            Scalar cl2 = polar.GetCl();
            cd_2d_(i) = polar.GetCd();
            dcl2d_dalpha_(i) = polar.GetDclDalphaPerRadian();
            dcd2d_dalpha_(i) = polar.GetDcdDalphaPerRadian();
            r(i) = cl_vlm(i) - cl2;
        }
    }
    return r;
}

// =============================================================================
// Forward chain: propagate a perturbation through VLM + Scirocco linearization
// =============================================================================

ViscousVlmCoupler::ForwardResponse ViscousVlmCoupler::RunForwardChain(
    const VecX& dn, const Vec3& dv_inf, const VecX& dtwist_explicit) {
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_vort = mesh.GetNumVortices();
    const Index n_strips = mesh.GetNumStrips();
    const MatX& vim = solver_.GetVelocityInfluenceMatrix();

    Vec3 v_inf_base = VortexLatticeSolver::ComputeFreestreamDirection(flight_);

    // -----------------------------------------------------------------
    // δR (VLM residual) = δAIC·Γ − δb, with:
    //   δAIC[i,j] = VIM[i,j] · δn_i        (VIM doesn't depend on twist)
    //   δb[i]     = −(v_inf + δv_inf) · n_i  →  δb[i] = −δv_inf · n_i
    // Hence  δR_i = δn_i · v_total_i  +  v_inf direction-sensitivity term:
    //
    //   δR_i =  δn_i · (v_inf + v_ind_i)
    //         + n_i  · δv_inf
    //
    // (The first term is the AIC perturbation collapsing to v_total at the
    // control point; the second is the RHS shift from the freestream tilt.)
    // -----------------------------------------------------------------
    VecX dR(n_vort);
    for (Index i = 0; i < n_vort; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);
        Vec3 v_total(v_inf_base.x() + last_solution_.induced_velocity(i, 0),
                     v_inf_base.y() + last_solution_.induced_velocity(i, 1),
                     v_inf_base.z() + last_solution_.induced_velocity(i, 2));
        Vec3 dni(dn(3 * i + 0), dn(3 * i + 1), dn(3 * i + 2));
        dR(i) = dni.dot(v_total) + vor.normal_at_control_point.dot(dv_inf);
    }

    // δΓ = AIC^{-1} · (−δR)   via cached LU factorization.
    VecX dGamma = solver_.SolveForward(-dR);

    // δv_ind_v = VIM[v, :] · δΓ  (3-component per vortex).
    // Pre-compute as a (3*N_vort) vector for ease.
    VecX dv_ind(3 * n_vort);
    for (Index i = 0; i < n_vort; ++i) {
        Scalar vx = 0.0;
        Scalar vy = 0.0;
        Scalar vz = 0.0;
        for (Index j = 0; j < n_vort; ++j) {
            Scalar gj = dGamma(j);
            vx += vim(3 * i + 0, j) * gj;
            vy += vim(3 * i + 1, j) * gj;
            vz += vim(3 * i + 2, j) * gj;
        }
        dv_ind(3 * i + 0) = vx;
        dv_ind(3 * i + 1) = vy;
        dv_ind(3 * i + 2) = vz;
    }

    // -----------------------------------------------------------------
    // δforce_v = 2*δΓ_v*(v_eff_v × dl_v) + 2*Γ_v*(δv_eff_v × dl_v)
    // δv_eff_v = δv_inf + δv_ind_v   (note: δv_inf is global; δv_ind_v is
    // per-vortex from above)
    // δdl_v = 0  (twist doesn't move bound vortex endpoints)
    //
    // From per-vortex δforces, aggregate per strip → δcf_strip; project
    // through BodyToStability(α∞) to get δcl_strip. Also aggregate the
    // whole-aircraft integral for δCL_total / δCD_total.
    // -----------------------------------------------------------------
    ForwardResponse out;
    out.dCl_strip = VecX::Zero(n_strips);
    out.dAlphaEff = VecX::Zero(n_strips);
    out.dCL_total = 0.0;
    out.dCD_total = 0.0;

    const ReferenceQuantities& ref = solver_.GetReferenceQuantities();
    const Scalar S_ref = std::max(ref.reference_area, Scalar(1e-30));
    const Scalar ca = std::cos(flight_.alpha);
    const Scalar sa = std::sin(flight_.alpha);

    Vec3 total_dforce = Vec3::Zero();

    // For α∞ perturbations the body→stability rotation also rotates the
    // base per-strip cf, contributing an extra term to δCl_strip. dα∞ is
    // encoded in dv_inf: with β=0, dv_inf = (-sin α∞, 0, cos α∞)·dα∞.
    Scalar dalpha_inf =
        (std::abs(ca) > Scalar(1e-12)) ? (dv_inf.z() / ca) : Scalar(0.0);

    for (Index strip_idx = 0; strip_idx < n_strips; ++strip_idx) {
        const Strip& strip = mesh.GetStrip(strip_idx);
        Vec3 strip_dforce = Vec3::Zero();

        for (Index k = 0; k < strip.num_vortices; ++k) {
            Index v = strip.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);

            Vec3 v_ind(last_solution_.induced_velocity(v, 0),
                       last_solution_.induced_velocity(v, 1),
                       last_solution_.induced_velocity(v, 2));
            Vec3 v_eff = v_inf_base + v_ind;
            Vec3 dl = vor.bound_vortex_right - vor.bound_vortex_left;

            Vec3 dv_eff(dv_inf.x() + dv_ind(3 * v + 0),
                        dv_inf.y() + dv_ind(3 * v + 1),
                        dv_inf.z() + dv_ind(3 * v + 2));

            Scalar gamma_v = last_solution_.circulation(v);
            Scalar dgamma_v = dGamma(v);

            Vec3 df = 2.0 * dgamma_v * v_eff.cross(dl) +
                      2.0 * gamma_v * dv_eff.cross(dl);
            strip_dforce += df;
        }

        // Per-strip δCl: normalize by strip area, project to stability axis.
        Scalar strip_area = strip.chord_center * strip.width;
        if (strip_area > Scalar(1e-30)) {
            Scalar dcf_x = strip_dforce.x() / strip_area;
            Scalar dcf_z = strip_dforce.z() / strip_area;
            out.dCl_strip(strip_idx) = dcf_z * ca - dcf_x * sa;

            // Body→stability rotation derivative on the base cf when α∞
            // changes. Without this term, ∂R/∂α∞ is missing a δCl_VLM
            // contribution and the implicit-function chain is biased.
            const auto& surf_loads =
                last_results_.surface_loads[strip.surface_index];
            for (const auto& sl : surf_loads.strip_loads) {
                if (sl.strip_index == strip_idx) {
                    out.dCl_strip(strip_idx) +=
                        dalpha_inf * (-sl.cf.z() * sa - sl.cf.x() * ca);
                    break;
                }
            }
        }

        // α_eff at the strip's LE vortex.
        Index iv = strip.first_vortex_index;
        Scalar vx_total =
            v_inf_base.x() + last_solution_.induced_velocity(iv, 0);
        Scalar vz_total =
            v_inf_base.z() + last_solution_.induced_velocity(iv, 2);
        Scalar denom = vx_total * vx_total + vz_total * vz_total;
        if (denom > Scalar(1e-30)) {
            Scalar dvx = dv_inf.x() + dv_ind(3 * iv + 0);
            Scalar dvz = dv_inf.z() + dv_ind(3 * iv + 2);
            // δ(atan2(z, x)) = (x*δz − z*δx) / (x² + z²)
            Scalar dalpha_body = (vx_total * dvz - vz_total * dvx) / denom;
            // α_eff = α_body − twist  ⇒  δα_eff = δα_body − δtwist
            out.dAlphaEff(strip_idx) = dalpha_body - dtwist_explicit(strip_idx);
        }

        total_dforce += strip_dforce;
    }

    // Integrated δCL, δCDi via BodyToStability at α∞. The rotation itself
    // depends on α∞ only when δv_inf is nonzero (α∞ perturbation); since
    // the strip-wise force aggregation already captured the body-axis
    // sensitivity, the rotation derivative is handled here as the final
    // step. For pure-Δα perturbations (dv_inf = 0) the rotation
    // contribution vanishes because cf_x and cf_z don't carry an α-rotation
    // term.
    Scalar dCF_x = total_dforce.x() / S_ref;
    Scalar dCF_z = total_dforce.z() / S_ref;
    out.dCL_total = dCF_z * ca - dCF_x * sa;
    Scalar dCDi_total = dCF_x * ca + dCF_z * sa;

    // For α∞ perturbations the body-to-stability rotation also rotates the
    // base integrated CF → extra term:
    //   d(CL)/dα = ... + (-CF_z·sin α - CF_x·cos α)·dα∞
    //   d(CD)/dα = ... + (-CF_x·sin α + CF_z·cos α)·dα∞
    // (dalpha_inf was computed above for the same reason in per-strip Cl.)
    Scalar CF_x = last_results_.CF.x();
    Scalar CF_z = last_results_.CF.z();
    out.dCL_total += dalpha_inf * (-CF_z * sa - CF_x * ca);
    dCDi_total += dalpha_inf * (-CF_x * sa + CF_z * ca);

    // Viscous drag contribution: δCDv = sum_i (dCd_2D/dα * δα_eff * area_i)
    // / S_ref. dCd_2D/dα is cached from the most recent Scirocco solve at
    // each strip.
    Scalar dCDv_total = 0.0;
    for (Index strip_idx = 0; strip_idx < n_strips; ++strip_idx) {
        const Strip& strip = mesh.GetStrip(strip_idx);
        Scalar strip_area = strip.chord_center * strip.width;
        dCDv_total +=
            dcd2d_dalpha_(strip_idx) * out.dAlphaEff(strip_idx) * strip_area;
    }
    dCDv_total /= S_ref;

    out.dCD_total = dCDi_total + dCDv_total;
    return out;
}

MatX ViscousVlmCoupler::BuildJacobianAnalytical() {
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_vort = mesh.GetNumVortices();
    const Index n_strips = mesh.GetNumStrips();
    MatX J(n_strips, n_strips);

    VecX dn = VecX::Zero(3 * n_vort);
    VecX dtwist = VecX::Zero(n_strips);
    Vec3 dv_inf_zero = Vec3::Zero();

    for (Index j = 0; j < n_strips; ++j) {
        // Build the normal perturbation for a unit Δα at strip j.
        const Strip& strip_j = mesh.GetStrip(j);
        dn.setZero();
        for (Index k = 0; k < strip_j.num_vortices; ++k) {
            Index v = strip_j.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);
            Vec3 dn_v = NormalDerivativeWrtTwist(vor, strip_j.twist_angle,
                                                 strip_j.trefftz_normal_y,
                                                 strip_j.trefftz_normal_z);
            dn.segment<3>(3 * v) = dn_v;
        }
        // The explicit α_eff contribution: changing Δα_j directly subtracts
        // 1 from α_eff at strip j (since α_eff = α_body − twist).
        dtwist.setZero();
        dtwist(j) = 1.0;

        ForwardResponse resp = RunForwardChain(dn, dv_inf_zero, dtwist);

        // J_ij = δR_i/δΔα_j = δCl_VLM,i − (dCl_2D,i/dα) · δα_eff,i
        for (Index i = 0; i < n_strips; ++i) {
            J(i, j) = resp.dCl_strip(i) - dcl2d_dalpha_(i) * resp.dAlphaEff(i);
        }
    }
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

        // Analytical Jacobian — one forward chain per unknown, reusing the
        // cached AIC LU. No additional VLM or Scirocco solves.
        MatX J = BuildJacobianAnalytical();

        Eigen::FullPivLU<MatX> lu(J);
        VecX du = lu.solve(-r);

        // Step-size limit: avoid overshooting into non-physical fixed
        // points when J is poorly conditioned. Scale uniformly so that
        // ‖Δu‖∞ ≤ options_.max_step_rad.
        Scalar max_abs = du.lpNorm<Eigen::Infinity>();
        if (max_abs > options_.max_step_rad) {
            du *= options_.max_step_rad / max_abs;
        }
        delta_alpha_ += options_.damping * du;
    }

    // Final evaluation at the converged Δα to refresh results.
    EvaluateVlm(cl_vlm, alpha_eff);
    (void)ComputeResidual(alpha_eff, cl_vlm);

    // Viscous drag bookkeeping: CD_viscous = Σ (Cd_2D,i * area_i) / S_ref.
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
    last_results_.CD = last_results_.CDi_nearfield + cdv;

    out.results = last_results_;
    out.delta_alpha = delta_alpha_;
    return out;
}

// =============================================================================
// Coupled total derivatives — fully analytical
// =============================================================================

CoupledGradients ViscousVlmCoupler::ComputeAlphaDerivatives() {
    CoupledGradients g;
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_vort = mesh.GetNumVortices();
    const Index n_strips = mesh.GetNumStrips();

    // Refresh base state so caches reflect the converged Δα.
    VecX cl_base(n_strips), alpha_eff_base(n_strips);
    EvaluateVlm(cl_base, alpha_eff_base);
    (void)ComputeResidual(alpha_eff_base, cl_base);

    g.CL = last_results_.CL;
    g.CD = last_results_.CDi_nearfield + last_results_.CDv;

    // -----------------------------------------------------------------
    // Build J = ∂R/∂u and ∂F/∂u  (both via analytical forward chain).
    // -----------------------------------------------------------------
    MatX J(n_strips, n_strips);
    VecX dCL_du = VecX::Zero(n_strips);
    VecX dCD_du = VecX::Zero(n_strips);

    VecX dn = VecX::Zero(3 * n_vort);
    VecX dtwist = VecX::Zero(n_strips);

    for (Index j = 0; j < n_strips; ++j) {
        const Strip& strip_j = mesh.GetStrip(j);
        dn.setZero();
        for (Index k = 0; k < strip_j.num_vortices; ++k) {
            Index v = strip_j.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);
            Vec3 dn_v = NormalDerivativeWrtTwist(vor, strip_j.twist_angle,
                                                 strip_j.trefftz_normal_y,
                                                 strip_j.trefftz_normal_z);
            dn.segment<3>(3 * v) = dn_v;
        }
        dtwist.setZero();
        dtwist(j) = 1.0;

        ForwardResponse resp = RunForwardChain(dn, Vec3::Zero(), dtwist);
        for (Index i = 0; i < n_strips; ++i) {
            J(i, j) = resp.dCl_strip(i) - dcl2d_dalpha_(i) * resp.dAlphaEff(i);
        }
        dCL_du(j) = resp.dCL_total;
        dCD_du(j) = resp.dCD_total;
    }

    // -----------------------------------------------------------------
    // ∂R/∂α∞ and ∂F/∂α∞ via a single forward chain with dv_inf = dα∞ direction.
    // -----------------------------------------------------------------
    double v_inf_arr[3];
    double dv_inf_arr[3];
    gradients::compute_freestream_direction(static_cast<double>(flight_.alpha),
                                            static_cast<double>(flight_.beta),
                                            v_inf_arr);
    gradients::compute_freestream_alpha_gradient(
        static_cast<double>(flight_.alpha), static_cast<double>(flight_.beta),
        dv_inf_arr);
    Vec3 dv_inf_alpha(static_cast<Scalar>(dv_inf_arr[0]),
                      static_cast<Scalar>(dv_inf_arr[1]),
                      static_cast<Scalar>(dv_inf_arr[2]));

    // α∞ perturbation has no explicit twist contribution (twist is fixed).
    VecX zero_twist = VecX::Zero(n_strips);
    VecX zero_dn = VecX::Zero(3 * n_vort);
    ForwardResponse resp_alpha =
        RunForwardChain(zero_dn, dv_inf_alpha, zero_twist);

    VecX dR_dalpha(n_strips);
    for (Index i = 0; i < n_strips; ++i) {
        dR_dalpha(i) = resp_alpha.dCl_strip(i) -
                       dcl2d_dalpha_(i) * resp_alpha.dAlphaEff(i);
    }
    Scalar dCL_dalpha_partial = resp_alpha.dCL_total;
    Scalar dCD_dalpha_partial = resp_alpha.dCD_total;

    // -----------------------------------------------------------------
    // Forward sensitivity:  J · (du/dα) = −∂R/∂α
    // Total derivative:     dF/dα = ∂F/∂α + ∂F/∂u · (du/dα)
    // -----------------------------------------------------------------
    Eigen::FullPivLU<MatX> lu(J);
    VecX du_dalpha = lu.solve(-dR_dalpha);

    g.dCL_dalpha = dCL_dalpha_partial + dCL_du.dot(du_dalpha);
    g.dCD_dalpha = dCD_dalpha_partial + dCD_du.dot(du_dalpha);

    return g;
}

}  // namespace pegasus::viscous
