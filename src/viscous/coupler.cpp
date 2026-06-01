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
#include <stdexcept>

#include <Eigen/Dense>

namespace pegasus::viscous {

namespace {

/// Local-frame triad at a strip. Built from strip geometry only — no flow
/// state — so it captures sweep (via spanwise_unit_vector) and dihedral
/// (via trefftz_normal components) consistently with how Pegasus places
/// panel normals.
struct StripLocalFrame {
    Vec3 c_hat;  ///< chord direction (LE → TE)
    Vec3 s_hat;  ///< spanwise direction (along bound vortex)
    Vec3 n_hat;  ///< lift normal (c_hat × s_hat, normalized)
};

/// d/dtwist of the normalized control-point normal at a vortex, evaluated
/// at the current strip twist. Same closed form used by v1 — the chord
/// direction is the only twist-dependent term in
/// normal = normalize(chord_dir × bound_dir).
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

StripLocalFrame ComputeStripLocalFrame(const Strip& strip) {
    Scalar t = strip.twist_angle;
    Scalar ny = strip.trefftz_normal_y;
    Scalar nz = strip.trefftz_normal_z;
    StripLocalFrame f;
    f.c_hat = Vec3(std::cos(t), -std::sin(t) * ny, -std::sin(t) * nz);
    f.s_hat = strip.spanwise_unit_vector;
    Vec3 cross = f.c_hat.cross(f.s_hat);
    Scalar mag = cross.norm();
    if (mag > Scalar(1e-30)) {
        f.n_hat = cross / mag;
    } else {
        f.n_hat = Vec3(0, 0, 1);
    }
    return f;
}

/// Per-strip effective angle of attack used as the Scirocco α input.
///
/// Compose two pieces:
///   1. α_body_local: project v_total onto the strip's chord-up plane
///      (drop the spanwise component), take atan2 in that plane. This
///      correctly handles sweep and dihedral via the strip's 3D triad.
///   2. Subtract the strip's effective twist (after the spline-injected
///      Δα offset). The minus sign makes the Newton Jacobian of the
///      Cl-mismatch residual monotonic in Δα with the physically-
///      intended fixed point near Δα = 0; the geometric "+twist" sign
///      lands at the spurious lift-free root.
///
/// For a planar wing this collapses to `α_body − twist`, matching the
/// original (working) v1 convention exactly.
Scalar ComputeStripAlphaEffLocal(const Strip& strip,
                                 const StripLocalFrame& frame,
                                 const Vec3& v_total) {
    Vec3 v_proj = v_total - (v_total.dot(frame.s_hat)) * frame.s_hat;
    Scalar v_c = v_proj.dot(frame.c_hat);
    Scalar v_n = v_proj.dot(frame.n_hat);
    Scalar alpha_body_local = Scalar(0.0);
    if (!(std::abs(v_c) < Scalar(1e-30) && std::abs(v_n) < Scalar(1e-30))) {
        alpha_body_local = std::atan2(v_n, v_c);
    }
    // The atan2 projection above gives α_body + twist for a planar wing,
    // so subtracting 2·twist yields the v1 "α_body − twist" sign. For
    // swept / dihedral wings this still puts the residual on the
    // monotonic branch that has the physical root at Δα ≈ 0.
    return alpha_body_local - Scalar(2.0) * strip.twist_angle;
}

}  // namespace

ViscousVlmCoupler::ViscousVlmCoupler(VortexLatticeSolver& solver,
                                     const FlightCondition& flight,
                                     const CouplerOptions& options)
    : solver_(solver), flight_(flight), options_(options) {
    if (options_.stations.size() < 2) {
        throw std::invalid_argument(
            "ViscousVlmCoupler: at least 2 coupling stations are required.");
    }
    const VortexLattice& mesh = solver_.GetMesh();
    interp_ = BuildStationInterpolation(mesh, options_.stations);

    const Index n_strips = mesh.GetNumStrips();
    const Index n_stations = static_cast<Index>(options_.stations.size());

    baseline_twist_.resize(n_strips);
    for (Index i = 0; i < n_strips; ++i) {
        baseline_twist_(i) = mesh.GetStrip(i).twist_angle;
    }
    delta_alpha_ = VecX::Zero(n_stations);
    cd_2d_station_ = VecX::Zero(n_stations);
    dcl2d_dalpha_station_ = VecX::Zero(n_stations);
    dcd2d_dalpha_station_ = VecX::Zero(n_stations);

    station_polars_.reserve(static_cast<std::size_t>(n_stations));
    for (const CouplingStation& s : options_.stations) {
        station_polars_.emplace_back(s.coords, s.num_panels);
    }
}

void ViscousVlmCoupler::ApplyTwistsAndRefresh() {
    // Spread per-station Δα to per-strip twist offsets via the spline op.
    VecX twist_offset = interp_.to_strip * delta_alpha_;
    VortexLattice& mesh = solver_.GetMeshMutable();
    const Index n = mesh.GetNumStrips();
    for (Index i = 0; i < n; ++i) {
        mesh.GetStrip(i).twist_angle = baseline_twist_(i) + twist_offset(i);
    }
    mesh.RecalculateNormals();
    solver_.RebuildAIC();
}

void ViscousVlmCoupler::EvaluateVlm(VecX& cl_local_strip,
                                    VecX& alpha_eff_strip) {
    ApplyTwistsAndRefresh();
    last_solution_ = solver_.Solve(flight_);
    last_results_ =
        force_calc_.ComputeForces(solver_.GetReferenceQuantities(),
                                  solver_.GetMesh(), last_solution_, flight_);

    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_strips = mesh.GetNumStrips();
    cl_local_strip.resize(n_strips);
    alpha_eff_strip.resize(n_strips);

    Vec3 v_inf = VortexLatticeSolver::ComputeFreestreamDirection(flight_);

    // For each strip: total body-axis force from its vortices; then
    //   Cl_VLM = stability-axis lift (cf.z·cos α∞ − cf.x·sin α∞), matching
    //            Pegasus's per-strip Cl convention exactly so the residual
    //            is well-defined and has the physical fixed point.
    //   α_eff = local-frame projection (chord-up plane) of v_total at the
    //           strip's LE vortex, then subtract 2·twist (see
    //           ComputeStripAlphaEffLocal for sign rationale).
    const Scalar ca = std::cos(flight_.alpha);
    const Scalar sa = std::sin(flight_.alpha);
    for (Index strip_idx = 0; strip_idx < n_strips; ++strip_idx) {
        const Strip& strip = mesh.GetStrip(strip_idx);
        StripLocalFrame frame = ComputeStripLocalFrame(strip);

        Vec3 total_force = Vec3::Zero();
        for (Index k = 0; k < strip.num_vortices; ++k) {
            Index v = strip.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);
            Vec3 v_ind(last_solution_.induced_velocity(v, 0),
                       last_solution_.induced_velocity(v, 1),
                       last_solution_.induced_velocity(v, 2));
            Vec3 v_eff = v_inf + v_ind;
            Vec3 dl = vor.bound_vortex_right - vor.bound_vortex_left;
            Scalar gamma = last_solution_.circulation(v);
            total_force += 2.0 * gamma * v_eff.cross(dl);
        }

        Scalar strip_area = strip.chord_center * strip.width;
        if (strip_area > Scalar(1e-30)) {
            Scalar cf_x = total_force.x() / strip_area;
            Scalar cf_z = total_force.z() / strip_area;
            cl_local_strip(strip_idx) = cf_z * ca - cf_x * sa;
        } else {
            cl_local_strip(strip_idx) = 0.0;
        }

        Index iv = strip.first_vortex_index;
        Vec3 v_ind_iv(last_solution_.induced_velocity(iv, 0),
                      last_solution_.induced_velocity(iv, 1),
                      last_solution_.induced_velocity(iv, 2));
        Vec3 v_total = v_inf + v_ind_iv;
        alpha_eff_strip(strip_idx) =
            ComputeStripAlphaEffLocal(strip, frame, v_total);
    }
}

VecX ViscousVlmCoupler::ComputeResidual(const VecX& alpha_eff_strip,
                                        const VecX& cl_local_strip) {
    // Interpolate per-strip α_eff and Cl_local onto stations.
    VecX alpha_eff_station = interp_.to_station * alpha_eff_strip;
    VecX cl_local_station = interp_.to_station * cl_local_strip;

    // Cl-mismatch residual: R_k = Cl_VLM_stability_station − Cl_2D(α_eff_k).
    // With α_eff carrying the −twist sign convention (see
    // ComputeStripAlphaEffLocal), the residual is monotonic in Δα with the
    // single fixed point near the physical (small) Δα solution. Scirocco
    // is driven in fixed-α mode at each station's interpolated α_eff;
    // its fixed-Cl mode is unreliable on this code path.
    const Index n_stations = static_cast<Index>(station_polars_.size());
    VecX r(n_stations);
    cd_2d_station_.resize(n_stations);
    dcl2d_dalpha_station_.resize(n_stations);
    dcd2d_dalpha_station_.resize(n_stations);

    for (Index k = 0; k < n_stations; ++k) {
        auto& polar = station_polars_[static_cast<std::size_t>(k)];
        Scalar re = options_.stations[static_cast<std::size_t>(k)].reynolds;
        bool ok =
            polar.SolveFixedAlpha(alpha_eff_station(k), re, options_.mach);
        if (!ok) {
            Scalar cl2 = 2.0 * kPi * alpha_eff_station(k);
            r(k) = cl_local_station(k) - cl2;
            cd_2d_station_(k) = 0.0;
            dcl2d_dalpha_station_(k) = 2.0 * kPi;
            dcd2d_dalpha_station_(k) = 0.0;
        } else {
            r(k) = cl_local_station(k) - polar.GetCl();
            cd_2d_station_(k) = polar.GetCd();
            dcl2d_dalpha_station_(k) = polar.GetDclDalphaPerRadian();
            dcd2d_dalpha_station_(k) = polar.GetDcdDalphaPerRadian();
        }
    }
    return r;
}

// =============================================================================
// Analytical forward chain
// =============================================================================
//
// For a perturbation direction (dtwist_strip, dv_inf), propagate analytically
// through the VLM solve + local-frame projections + spline interp + Scirocco
// linearization to get the response of all coupled quantities:
//
//   dtwist_strip ─►  δn_v at each vortex (closed form, NormalDeriv)
//                   δR_VLM_v = δn_v · v_total_v  +  n_v · δv_inf
//                   δΓ via cached AIC LU
//                   δv_ind = VIM · δΓ
//                   δforce_v = 2·δΓ·(v_eff×dl) + 2·Γ·(δv_eff×dl)
//                   δstrip_force = Σ_v δforce_v
//                   δc_hat, δn_hat from δtwist_strip (closed form)
//                   δCl_local_strip = (δforce·n_hat + force·δn_hat)/area
//                   δα_eff_strip via atan2 chain rule (incl. δv,δn̂,δĉ)
//                                 + −2·δtwist_strip
//                   to_station ·  ⟶  per-station δCl, δα_eff
//                   Σ δforce, BodyToStability ⟶ δCL, δCDi
//                   (dCd_2D/dα)·δα_eff_station, to_strip, sum ⟶ δCDv
//
// All linear ops are precomputed/cached at coupler state. The δn / δn_hat
// / δc_hat closed forms exploit that twist changes only the panel's
// rotation about its spanwise axis (bound vortex positions stay fixed).
ViscousVlmCoupler::ForwardResponse ViscousVlmCoupler::RunForwardChain(
    const VecX& dtwist_strip, const Vec3& dv_inf) {
    const VortexLattice& mesh = solver_.GetMesh();
    const MatX& vim = solver_.GetVelocityInfluenceMatrix();
    const Index n_vort = mesh.GetNumVortices();
    const Index n_strips = mesh.GetNumStrips();
    const Index n_stations = static_cast<Index>(station_polars_.size());

    Vec3 v_inf_base = VortexLatticeSolver::ComputeFreestreamDirection(flight_);

    // δn_cp (per vortex, 3 components flattened). Only the vortices in
    // strips with non-zero δtwist_strip contribute. The closed-form
    // NormalDerivativeWrtTwist gives the derivative of normal_at_cp w.r.t.
    // the strip's twist_angle at the current state.
    VecX dn_cp = VecX::Zero(3 * n_vort);
    for (Index strip_idx = 0; strip_idx < n_strips; ++strip_idx) {
        Scalar dt = dtwist_strip(strip_idx);
        if (std::abs(dt) < Scalar(1e-30)) {
            continue;
        }
        const Strip& strip = mesh.GetStrip(strip_idx);
        for (Index k = 0; k < strip.num_vortices; ++k) {
            Index v = strip.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);
            Vec3 dn_v = NormalDerivativeWrtTwist(vor, strip.twist_angle,
                                                 strip.trefftz_normal_y,
                                                 strip.trefftz_normal_z);
            dn_cp.segment<3>(3 * v) = dt * dn_v;
        }
    }

    // VLM residual perturbation: δR_i = δn_i · v_total_i + n_i · δv_inf.
    VecX dR(n_vort);
    for (Index i = 0; i < n_vort; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);
        Vec3 v_total(v_inf_base.x() + last_solution_.induced_velocity(i, 0),
                     v_inf_base.y() + last_solution_.induced_velocity(i, 1),
                     v_inf_base.z() + last_solution_.induced_velocity(i, 2));
        Vec3 dni(dn_cp(3 * i + 0), dn_cp(3 * i + 1), dn_cp(3 * i + 2));
        dR(i) = dni.dot(v_total) + vor.normal_at_control_point.dot(dv_inf);
    }

    // δΓ via cached AIC LU.
    VecX dGamma = solver_.SolveForward(-dR);

    // δv_ind = VIM · δΓ (precompute for all vortices).
    VecX dv_ind = VecX::Zero(3 * n_vort);
    for (Index i = 0; i < n_vort; ++i) {
        Scalar vx = 0.0, vy = 0.0, vz = 0.0;
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

    // Per-strip aggregation: total δforce_body, base strip force, base
    // local frame, then δc_hat/δn_hat from δtwist_strip and projection.
    ForwardResponse out;
    VecX dCl_local_strip = VecX::Zero(n_strips);
    VecX dAlphaEff_strip = VecX::Zero(n_strips);
    Vec3 total_dforce = Vec3::Zero();
    const ReferenceQuantities& ref = solver_.GetReferenceQuantities();
    const Scalar S_ref = std::max(ref.reference_area, Scalar(1e-30));
    const Scalar ca = std::cos(flight_.alpha);
    const Scalar sa = std::sin(flight_.alpha);
    Scalar dalpha_inf =
        (std::abs(ca) > Scalar(1e-12)) ? (dv_inf.z() / ca) : Scalar(0.0);

    for (Index si = 0; si < n_strips; ++si) {
        const Strip& strip = mesh.GetStrip(si);
        Scalar t = strip.twist_angle;
        Scalar ny = strip.trefftz_normal_y;
        Scalar nz = strip.trefftz_normal_z;
        Vec3 c_hat(std::cos(t), -std::sin(t) * ny, -std::sin(t) * nz);
        Vec3 s_hat = strip.spanwise_unit_vector;
        Vec3 cross_cs = c_hat.cross(s_hat);
        Scalar mag = cross_cs.norm();
        Vec3 n_hat = (mag > Scalar(1e-30)) ? (cross_cs / mag) : Vec3(0, 0, 1);

        Scalar dt = dtwist_strip(si);

        // δc_hat and δn_hat from δtwist_strip; δs_hat = 0 (spanwise
        // direction doesn't depend on twist).
        Vec3 d_c_hat =
            Vec3(-std::sin(t), -std::cos(t) * ny, -std::cos(t) * nz) * dt;
        Vec3 d_cross = d_c_hat.cross(s_hat);
        Vec3 d_n_hat = Vec3::Zero();
        if (mag > Scalar(1e-30)) {
            d_n_hat = (d_cross - n_hat * (n_hat.dot(d_cross))) / mag;
        }

        // Strip base force and its perturbation.
        Vec3 strip_force = Vec3::Zero();
        Vec3 strip_dforce = Vec3::Zero();
        for (Index k = 0; k < strip.num_vortices; ++k) {
            Index v = strip.first_vortex_index + k;
            const HorseshoeVortex& vor = mesh.GetVortex(v);
            Vec3 v_ind(last_solution_.induced_velocity(v, 0),
                       last_solution_.induced_velocity(v, 1),
                       last_solution_.induced_velocity(v, 2));
            Vec3 v_eff = v_inf_base + v_ind;
            Vec3 dl = vor.bound_vortex_right - vor.bound_vortex_left;
            Scalar gamma = last_solution_.circulation(v);
            Scalar dgamma = dGamma(v);
            Vec3 dv_eff(dv_inf.x() + dv_ind(3 * v + 0),
                        dv_inf.y() + dv_ind(3 * v + 1),
                        dv_inf.z() + dv_ind(3 * v + 2));
            strip_force += 2.0 * gamma * v_eff.cross(dl);
            strip_dforce +=
                2.0 * dgamma * v_eff.cross(dl) + 2.0 * gamma * dv_eff.cross(dl);
        }

        Scalar strip_area = strip.chord_center * strip.width;
        if (strip_area > Scalar(1e-30)) {
            // Stability-axis δCl: project δforce by α∞ rotation, plus the
            // body-to-stability rotation derivative on the base strip cf
            // when α∞ itself is perturbed (carried by dalpha_inf).
            Scalar dcf_x = strip_dforce.x() / strip_area;
            Scalar dcf_z = strip_dforce.z() / strip_area;
            Scalar cf_x_base = 0.0, cf_z_base = 0.0;
            const auto& surf_loads =
                last_results_.surface_loads[strip.surface_index];
            for (const auto& sl : surf_loads.strip_loads) {
                if (sl.strip_index == si) {
                    cf_x_base = sl.cf.x();
                    cf_z_base = sl.cf.z();
                    break;
                }
            }
            dCl_local_strip(si) =
                dcf_z * ca - dcf_x * sa +
                dalpha_inf * (-cf_z_base * sa - cf_x_base * ca);
        }

        // δα_eff_strip via atan2 chain rule.
        // Base v_total at the strip's LE vortex.
        Index iv = strip.first_vortex_index;
        Vec3 v_ind_iv(last_solution_.induced_velocity(iv, 0),
                      last_solution_.induced_velocity(iv, 1),
                      last_solution_.induced_velocity(iv, 2));
        Vec3 v_total = v_inf_base + v_ind_iv;
        Vec3 v_proj = v_total - (v_total.dot(s_hat)) * s_hat;
        Scalar v_c = v_proj.dot(c_hat);
        Scalar v_n = v_proj.dot(n_hat);

        Vec3 dv_total(dv_inf.x() + dv_ind(3 * iv + 0),
                      dv_inf.y() + dv_ind(3 * iv + 1),
                      dv_inf.z() + dv_ind(3 * iv + 2));
        // δs_hat = 0, so δv_proj = δv_total − (δv_total·s_hat)·s_hat.
        Vec3 dv_proj = dv_total - (dv_total.dot(s_hat)) * s_hat;

        Scalar d_v_c = dv_proj.dot(c_hat) + v_proj.dot(d_c_hat);
        Scalar d_v_n = dv_proj.dot(n_hat) + v_proj.dot(d_n_hat);

        Scalar denom = v_c * v_c + v_n * v_n;
        Scalar d_atan = (denom > Scalar(1e-30))
                            ? ((v_c * d_v_n - v_n * d_v_c) / denom)
                            : Scalar(0.0);
        // α_eff = atan2(v_n, v_c) − 2·twist  (see ComputeStripAlphaEffLocal).
        dAlphaEff_strip(si) = d_atan - Scalar(2.0) * dt;

        total_dforce += strip_dforce;
    }

    // Interpolate per-strip to per-station.
    out.dCl_station = interp_.to_station * dCl_local_strip;
    out.dAlphaEff_station = interp_.to_station * dAlphaEff_strip;

    // Integrated δCL / δCDi via body-to-stability rotation. For α∞
    // perturbations, also add the rotation derivative on the base CF.
    Scalar dCF_x = total_dforce.x() / S_ref;
    Scalar dCF_z = total_dforce.z() / S_ref;
    out.dCL_total = dCF_z * ca - dCF_x * sa;
    Scalar dCDi_total = dCF_x * ca + dCF_z * sa;
    Scalar CF_x = last_results_.CF.x();
    Scalar CF_z = last_results_.CF.z();
    out.dCL_total += dalpha_inf * (-CF_z * sa - CF_x * ca);
    dCDi_total += dalpha_inf * (-CF_x * sa + CF_z * ca);

    // δCDv contribution from each station's dCd_2D/dα · δα_eff, distributed
    // to strips by to_strip, then strip-area-weighted into the global ref.
    // VecX d_cd_station = (dcd2d_dalpha_station_ array, elementwise) ×
    // out.dAlphaEff_station.
    VecX d_cd_station(n_stations);
    for (Index k = 0; k < n_stations; ++k) {
        d_cd_station(k) = dcd2d_dalpha_station_(k) * out.dAlphaEff_station(k);
    }
    VecX d_cd_strip = interp_.to_strip * d_cd_station;
    Scalar dCDv_total = 0.0;
    for (Index i = 0; i < n_strips; ++i) {
        const Strip& strip = mesh.GetStrip(i);
        dCDv_total += d_cd_strip(i) * strip.chord_center * strip.width;
    }
    dCDv_total /= S_ref;

    out.dCD_total = dCDi_total + dCDv_total;
    return out;
}

MatX ViscousVlmCoupler::BuildJacobianAnalytical() {
    const Index n_stations = static_cast<Index>(station_polars_.size());
    MatX J(n_stations, n_stations);

    // R_i = Cl_VLM_i − Cl_2D(α_eff_i), so
    //     ∂R_i/∂u_k = ∂Cl_i/∂u_k − (dCl_2D/dα)_i · ∂α_eff_i/∂u_k.
    for (Index k = 0; k < n_stations; ++k) {
        VecX dtwist = interp_.to_strip.col(k);
        ForwardResponse resp = RunForwardChain(dtwist, Vec3::Zero());
        for (Index i = 0; i < n_stations; ++i) {
            J(i, k) = resp.dCl_station(i) -
                      dcl2d_dalpha_station_(i) * resp.dAlphaEff_station(i);
        }
    }
    return J;
}

CouplerResult ViscousVlmCoupler::Solve() {
    CouplerResult out;
    const Index n_stations = static_cast<Index>(station_polars_.size());

    VecX cl_local_strip, alpha_eff_strip;
    for (int iter = 0; iter < options_.max_iterations; ++iter) {
        EvaluateVlm(cl_local_strip, alpha_eff_strip);
        VecX r = ComputeResidual(alpha_eff_strip, cl_local_strip);
        Scalar res_norm = r.lpNorm<Eigen::Infinity>();
        out.iterations = iter + 1;
        out.final_residual = res_norm;
        if (res_norm < options_.tolerance) {
            out.converged = true;
            break;
        }

        MatX J = BuildJacobianAnalytical();
        Eigen::FullPivLU<MatX> lu(J);
        VecX du = lu.solve(-r);

        Scalar max_abs = du.lpNorm<Eigen::Infinity>();
        if (max_abs > options_.max_step_rad) {
            du *= options_.max_step_rad / max_abs;
        }

        // Backtracking line search: try the full step first, then halve
        // until ‖R‖ decreases. With the geometrically correct local-frame
        // α_eff convention the Newton J is well-defined but its diagonal
        // can be small in magnitude (3D lift slope close to 2π), so a full
        // step can overshoot into a non-physical lift-free fixed point.
        // Line search keeps Newton on the physical solution branch.
        const VecX delta_alpha_save = delta_alpha_;
        Scalar alpha_ls = options_.damping;
        const Scalar r_norm = res_norm;
        VecX cl_ls, alpha_ls_strip, r_ls;
        bool accepted = false;
        for (int bt = 0; bt < 10; ++bt) {
            delta_alpha_ = delta_alpha_save + alpha_ls * du;
            EvaluateVlm(cl_ls, alpha_ls_strip);
            r_ls = ComputeResidual(alpha_ls_strip, cl_ls);
            if (r_ls.lpNorm<Eigen::Infinity>() <= r_norm) {
                accepted = true;
                break;
            }
            alpha_ls *= Scalar(0.5);
        }
        if (!accepted) {
            // Step couldn't reduce ‖R‖ even at α = damping·2^−10. Take the
            // smallest tried step anyway and let Newton try again next
            // iter from a perturbed state.
            delta_alpha_ = delta_alpha_save + alpha_ls * du;
        }
    }

    // Refresh at the final Δα so caches reflect the returned state.
    EvaluateVlm(cl_local_strip, alpha_eff_strip);
    (void)ComputeResidual(alpha_eff_strip, cl_local_strip);

    // Viscous drag bookkeeping. CDv is integrated from station Cd via
    // interpolation back onto strips, then strip_area-weighted into the
    // global reference area.
    const VortexLattice& mesh = solver_.GetMesh();
    VecX cd_2d_strip = interp_.to_strip * cd_2d_station_;
    Scalar ref_area = std::max(last_results_.reference_area, Scalar(1e-30));
    Scalar cdv = 0.0;
    for (Index i = 0; i < mesh.GetNumStrips(); ++i) {
        const Strip& strip = mesh.GetStrip(i);
        cdv += cd_2d_strip(i) * strip.chord_center * strip.width;
    }
    cdv /= ref_area;
    last_results_.CDv = cdv;
    last_results_.CD = last_results_.CDi_nearfield + cdv;

    out.results = last_results_;
    out.delta_alpha_station = delta_alpha_;
    out.twist_offset_strip = interp_.to_strip * delta_alpha_;
    if (n_stations > 0 && !out.converged) {
        // Even if max_iterations was hit, still report the final state.
    }
    return out;
}

CoupledGradients ViscousVlmCoupler::ComputeAlphaDerivatives() {
    CoupledGradients g;
    const Index n_stations = static_cast<Index>(station_polars_.size());

    // Refresh base state so caches reflect the converged Δα.
    VecX cl_strip, alpha_strip;
    EvaluateVlm(cl_strip, alpha_strip);
    (void)ComputeResidual(alpha_strip, cl_strip);

    // Base integrated CL/CD (with viscous correction).
    const VortexLattice& mesh = solver_.GetMesh();
    Scalar ref_area = std::max(last_results_.reference_area, Scalar(1e-30));
    Scalar cdv_base = 0.0;
    {
        VecX cd_strip = interp_.to_strip * cd_2d_station_;
        for (Index i = 0; i < mesh.GetNumStrips(); ++i) {
            const Strip& s = mesh.GetStrip(i);
            cdv_base += cd_strip(i) * s.chord_center * s.width;
        }
        cdv_base /= ref_area;
    }
    g.CL = last_results_.CL;
    g.CD = last_results_.CDi_nearfield + cdv_base;

    // ────────────────────────────────────────────────────────────────────
    // Analytical: J = ∂R/∂u and ∂F/∂u via one forward chain per station.
    // ────────────────────────────────────────────────────────────────────
    MatX J(n_stations, n_stations);
    VecX dCL_du(n_stations);
    VecX dCD_du(n_stations);
    for (Index k = 0; k < n_stations; ++k) {
        VecX dtwist = interp_.to_strip.col(k);
        ForwardResponse resp = RunForwardChain(dtwist, Vec3::Zero());
        for (Index i = 0; i < n_stations; ++i) {
            J(i, k) = resp.dCl_station(i) -
                      dcl2d_dalpha_station_(i) * resp.dAlphaEff_station(i);
        }
        dCL_du(k) = resp.dCL_total;
        dCD_du(k) = resp.dCD_total;
    }

    // ────────────────────────────────────────────────────────────────────
    // Analytical: ∂R/∂α∞ and ∂F/∂α∞ via one forward chain with
    // dtwist_strip = 0 and dv_inf = df_inf/dα.
    // ────────────────────────────────────────────────────────────────────
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
    VecX zero_twist = VecX::Zero(mesh.GetNumStrips());
    ForwardResponse resp_alpha = RunForwardChain(zero_twist, dv_inf_alpha);
    VecX dR_dalpha(n_stations);
    for (Index k = 0; k < n_stations; ++k) {
        dR_dalpha(k) =
            resp_alpha.dCl_station(k) -
            dcl2d_dalpha_station_(k) * resp_alpha.dAlphaEff_station(k);
    }
    Scalar dCL_dalpha_partial = resp_alpha.dCL_total;
    Scalar dCD_dalpha_partial = resp_alpha.dCD_total;

    // Implicit function theorem: du/dα = −J^(−1)·∂R/∂α,  then
    // dF/dα = ∂F/∂α + ∂F/∂u · du/dα.
    Eigen::FullPivLU<MatX> lu(J);
    VecX du_dalpha = lu.solve(-dR_dalpha);
    g.dCL_dalpha = dCL_dalpha_partial + dCL_du.dot(du_dalpha);
    g.dCD_dalpha = dCD_dalpha_partial + dCD_du.dot(du_dalpha);
    return g;
}

// =============================================================================
// Coupled adjoint w.r.t. per-station airfoil coordinates
// =============================================================================
//
// For station-k coordinates x_k (2 × N matrix):
//     ∂F/∂x_k|_u : CL = 0;  CD = w_k · ∂Cd_2D_k/∂x_k
//                  with w_k = Σ_strip (to_strip)[strip, k] · area / S_ref
//     ∂R_j/∂x_k|_u : zero except for row j == k, where it equals
//                     −∂Cl_2D_k/∂x_k.
//
// Solve adjoint   J^T · λ_F = (∂F/∂u)^T   once per QoI:
//     dF/dx_k = ∂F/∂x_k|_u + λ_F[k] · ∂R_k/∂x_k
//             = ∂F/∂x_k|_u − λ_F[k] · ∂Cl_2D_k/∂x_k
//
// ∂F/∂u comes from the same RunForwardChain sweep used by
// ComputeAlphaDerivatives; ∂Cl_2D_k/∂x_k and ∂Cd_2D_k/∂x_k come from
// Scirocco's adjoint (GetDclDx / GetDcdDx, cached per station from the
// most recent SolveFixedAlpha).
CoupledShapeGradients ViscousVlmCoupler::ComputeShapeDerivatives() {
    CoupledShapeGradients g;
    const Index n_stations = static_cast<Index>(station_polars_.size());
    const VortexLattice& mesh = solver_.GetMesh();
    const Index n_strips = mesh.GetNumStrips();
    const Scalar S_ref = std::max(
        solver_.GetReferenceQuantities().reference_area, Scalar(1e-30));

    // Refresh base state and ensure Scirocco caches are populated.
    VecX cl_strip, alpha_strip;
    EvaluateVlm(cl_strip, alpha_strip);
    (void)ComputeResidual(alpha_strip, cl_strip);

    // Rebuild J and ∂F/∂u using the same forward sweep as
    // ComputeAlphaDerivatives.
    MatX J(n_stations, n_stations);
    VecX dCL_du(n_stations);
    VecX dCD_du(n_stations);
    for (Index k = 0; k < n_stations; ++k) {
        VecX dtwist = interp_.to_strip.col(k);
        ForwardResponse resp = RunForwardChain(dtwist, Vec3::Zero());
        for (Index i = 0; i < n_stations; ++i) {
            J(i, k) = resp.dCl_station(i) -
                      dcl2d_dalpha_station_(i) * resp.dAlphaEff_station(i);
        }
        dCL_du(k) = resp.dCL_total;
        dCD_du(k) = resp.dCD_total;
    }

    // Adjoint solves: J^T · λ_F = (∂F/∂u)^T, one per QoI.
    Eigen::FullPivLU<MatX> lu_t(J.transpose());
    VecX lambda_CL = lu_t.solve(dCL_du);
    VecX lambda_CD = lu_t.solve(dCD_du);

    // Per-station "viscous-drag weight" w_k = Σ_strip (to_strip)[strip,k] *
    // strip_area / S_ref. This is how strongly station k's drag enters the
    // integrated CDv.
    VecX w(n_stations);
    for (Index k = 0; k < n_stations; ++k) {
        Scalar sum = 0.0;
        for (Index i = 0; i < n_strips; ++i) {
            const Strip& s = mesh.GetStrip(i);
            sum += interp_.to_strip(i, k) * s.chord_center * s.width;
        }
        w(k) = sum / S_ref;
    }

    g.dCL_dcoords.resize(static_cast<std::size_t>(n_stations));
    g.dCD_dcoords.resize(static_cast<std::size_t>(n_stations));
    for (Index k = 0; k < n_stations; ++k) {
        auto& polar = station_polars_[static_cast<std::size_t>(k)];
        const Eigen::MatrixXd& dcl_dx = polar.GetDclDx();
        const Eigen::MatrixXd& dcd_dx = polar.GetDcdDx();

        // dF/dx = ∂F/∂x|_u − λ_F^T · ∂R/∂x. For x_k station-k coords:
        //   ∂R_j/∂x_k = −∂Cl_2D_k/∂x_k · δ_jk  (only row k, with a minus)
        // ⇒ the contribution from λ flips sign and becomes +λ[k]·dCl/dx.
        //
        // CL: ∂CL/∂x_k|_u = 0  (no direct effect on the integrated CL)
        //     dCL/dx_k = +λ_CL[k] · dCl_2D_k/dx_k
        Eigen::MatrixXd dCL_dx_k = lambda_CL(k) * dcl_dx;

        // CD: ∂CD/∂x_k|_u = w_k · dCd_2D_k/dx_k (the viscous-drag direct
        //     effect) and the implicit piece is +λ_CD[k]·dCl/dx_k.
        Eigen::MatrixXd dCD_dx_k = w(k) * dcd_dx + lambda_CD(k) * dcl_dx;

        g.dCL_dcoords[static_cast<std::size_t>(k)] = std::move(dCL_dx_k);
        g.dCD_dcoords[static_cast<std::size_t>(k)] = std::move(dCD_dx_k);
    }

    return g;
}

}  // namespace pegasus::viscous
