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

#include "vlm/ad/gradient_calculator.h"

#include "vlm/ad/adjoint_solver.h"
#include "vlm/ad/analytical_aic.h"
#include "vlm/ad/analytical_forces.h"
#include "vlm/ad/mesh_node_sensitivity.h"
#include "vlm/ad/mesh_sensitivity.h"
#include "vlm/aero/force_calculator.h"
#include "vlm/geometry/mesh_input.h"
#include "vlm/geometry/surface.h"

#include <cmath>
#include <vector>

namespace vlm::ad {

// ============================================================================
// Implementation Class
// ============================================================================

class GradientCalculator::Impl {
public:
    Aircraft aircraft_{"GradientCalc"};
    VortexLatticeSolver solver_;
    ForceCalculator force_calc_;
    bool aircraft_set_{false};
    bool mesh_set_{false};
    bool use_adjoint_{true};
    MeshInput mesh_input_;
    ReferenceQuantities ref_quantities_;

    // Cached solution for adjoint computations
    SolutionState cached_solution_;
    FlightCondition cached_condition_;
    bool solution_cached_{false};

    /// Solve and cache the solution at a flight condition
    void SolveAndCache(const FlightCondition& condition) {
        cached_solution_ = solver_.Solve(condition);
        cached_condition_ = condition;
        solution_cached_ = true;
    }

    /// Compute aerodynamic results at a flight condition (full solve)
    AerodynamicResults ComputeResults(const FlightCondition& condition) {
        SolutionState solution = solver_.Solve(condition);
        return force_calc_.ComputeForces(aircraft_, solver_.GetMesh(), solution,
                                         condition);
    }

    /// Compute RHS sensitivity ∂b/∂(parameter)
    /// The RHS is b_i = -n_i · V_eff where V_eff depends on flight condition
    VecX ComputeRHSSensitivity(const FlightCondition& condition,
                               int param_index) {
        const VortexLattice& mesh = solver_.GetMesh();
        const Index n = cached_solution_.circulation.size();
        VecX db_dparam(n);

        // Compute freestream direction gradient analytically
        Vec3 dV_dparam = Vec3::Zero();

        switch (param_index) {
            case 0: {  // alpha - analytical gradient
                double dv_dalpha[3];
                compute_freestream_alpha_gradient(condition.alpha,
                                                  condition.beta, dv_dalpha);
                dV_dparam = Vec3(dv_dalpha[0], dv_dalpha[1], dv_dalpha[2]);
                break;
            }
            case 1: {  // beta - analytical gradient
                double dv_dbeta[3];
                compute_freestream_beta_gradient(condition.alpha,
                                                 condition.beta, dv_dbeta);
                dV_dparam = Vec3(dv_dbeta[0], dv_dbeta[1], dv_dbeta[2]);
                break;
            }
            case 2:  // roll rate - freestream direction doesn't depend on p
            case 3:  // pitch rate - freestream direction doesn't depend on q
            case 4:  // yaw rate - freestream direction doesn't depend on r
            default:
                // dV_dparam stays zero; the rotation contribution is handled
                // below
                break;
        }

        // RHS sensitivity: ∂b_i/∂param = -n_i · (∂V/∂param)
        for (Index i = 0; i < n; ++i) {
            Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
            db_dparam(i) = -normal.dot(dV_dparam);
        }

        // For rotation rates, there's also a contribution from the rotation
        // velocity The RHS includes -n · (ω × r) where r is the control point
        // position
        if (param_index >= 2 && param_index <= 4) {
            const Vec3& ref = solver_.GetReferencePoint();

            for (Index i = 0; i < n; ++i) {
                Vec3 normal = mesh.GetVortex(i).normal_at_control_point;
                Vec3 cp = mesh.GetVortex(i).control_point;
                Vec3 r = cp - ref;

                Vec3 omega_dir;
                switch (param_index) {
                    case 2:
                        omega_dir = Vec3(1, 0, 0);
                        break;  // p (roll)
                    case 3:
                        omega_dir = Vec3(0, 1, 0);
                        break;  // q (pitch)
                    case 4:
                        omega_dir = Vec3(0, 0, 1);
                        break;  // r (yaw)
                    default:
                        break;
                }

                // Rotation velocity: V_rot = ω × r
                // ∂V_rot/∂ω = (direction) × r
                Vec3 dVrot = omega_dir.cross(r);
                db_dparam(i) -= normal.dot(dVrot);
            }
        }

        return db_dparam;
    }

    /// Extract mesh data into raw arrays for analytical kernels
    void ExtractMeshData(std::vector<double>& bound_left,
                         std::vector<double>& bound_right,
                         std::vector<double>& control_points,
                         std::vector<double>& moment_arm,
                         std::vector<double>& induced_velocity) const {
        const VortexLattice& mesh = solver_.GetMesh();
        const Index n = mesh.GetNumVortices();
        bound_left.resize(3 * n);
        bound_right.resize(3 * n);
        control_points.resize(3 * n);
        moment_arm.resize(3 * n);
        induced_velocity.resize(3 * n);

        Vec3 moment_ref = aircraft_.GetMomentReferencePoint();

        for (Index i = 0; i < n; ++i) {
            const HorseshoeVortex& vor = mesh.GetVortex(i);

            bound_left[3 * i + 0] = vor.bound_vortex_left.x();
            bound_left[3 * i + 1] = vor.bound_vortex_left.y();
            bound_left[3 * i + 2] = vor.bound_vortex_left.z();

            bound_right[3 * i + 0] = vor.bound_vortex_right.x();
            bound_right[3 * i + 1] = vor.bound_vortex_right.y();
            bound_right[3 * i + 2] = vor.bound_vortex_right.z();

            control_points[3 * i + 0] = vor.control_point.x();
            control_points[3 * i + 1] = vor.control_point.y();
            control_points[3 * i + 2] = vor.control_point.z();

            moment_arm[3 * i + 0] = vor.control_point.x() - moment_ref.x();
            moment_arm[3 * i + 1] = vor.control_point.y() - moment_ref.y();
            moment_arm[3 * i + 2] = vor.control_point.z() - moment_ref.z();

            if (cached_solution_.induced_velocity.rows() > i) {
                induced_velocity[3 * i + 0] =
                    cached_solution_.induced_velocity(i, 0);
                induced_velocity[3 * i + 1] =
                    cached_solution_.induced_velocity(i, 1);
                induced_velocity[3 * i + 2] =
                    cached_solution_.induced_velocity(i, 2);
            } else {
                induced_velocity[3 * i + 0] = 0.0;
                induced_velocity[3 * i + 1] = 0.0;
                induced_velocity[3 * i + 2] = 0.0;
            }
        }
    }

    /// Compute explicit term ∂J/∂param at fixed circulation using analytical
    /// derivatives
    AeroGradients ComputeExplicitDerivatives(int param_index) {
        const Index n = cached_solution_.circulation.size();
        if (n == 0) {
            return AeroGradients{};
        }

        // Extract mesh data
        std::vector<double> bound_left, bound_right, control_points, moment_arm,
            induced_velocity;
        ExtractMeshData(bound_left, bound_right, control_points, moment_arm,
                        induced_velocity);

        // Copy circulation to contiguous array
        std::vector<double> circulation(n);
        for (Index i = 0; i < n; ++i) {
            circulation[i] = cached_solution_.circulation(i);
        }

        double reference_area = aircraft_.GetReferenceArea();
        double d_result[6];

        AeroGradients grads;

        if (param_index == 0) {
            // Alpha derivative
            analytical_force_alpha_gradient(
                cached_condition_.alpha, cached_condition_.beta,
                circulation.data(), induced_velocity.data(), bound_left.data(),
                bound_right.data(), control_points.data(), moment_arm.data(),
                static_cast<int>(n), reference_area, d_result);
        } else if (param_index == 1) {
            // Beta derivative
            analytical_force_beta_gradient(
                cached_condition_.alpha, cached_condition_.beta,
                circulation.data(), induced_velocity.data(), bound_left.data(),
                bound_right.data(), control_points.data(), moment_arm.data(),
                static_cast<int>(n), reference_area, d_result);
        } else {
            // Rotation rates (p, q, r): The force kernel only uses V_inf +
            // V_induced and does not include V_rotation. The rotation effect on
            // forces is entirely captured through the implicit term (how Γ
            // changes when the RHS changes with rotation rate). The explicit
            // force term is zero.
            for (double& out_idx : d_result) {
                out_idx = 0.0;
            }
        }

        grads.dCL = d_result[0];
        grads.dCD = d_result[1];
        grads.dCY = d_result[2];
        grads.dCl = d_result[3];
        grads.dCm = d_result[4];
        grads.dCn = d_result[5];

        // CDi_ff (Trefftz plane drag) doesn't have an explicit
        // alpha/beta/rotation term because it's computed from circulation and
        // geometry only
        grads.dCDi = 0.0;

        return grads;
    }

    /// Compute derivatives using adjoint method
    /// dJ/d(param) = ∂J/∂param|_Γ + λᵀ·(∂b/∂param)
    AeroGradients AdjointDerivatives(const FlightCondition& condition,
                                     int param_index) {
        // Step 1: Solve at base condition and cache
        if (!solution_cached_ || condition.alpha != cached_condition_.alpha ||
            condition.beta != cached_condition_.beta) {
            SolveAndCache(condition);
        }

        // Step 2: Compute force gradients ∂J/∂Γ for all outputs
        ForceGradients force_grads = ComputeForceGradients(
            solver_.GetMesh(), cached_solution_, condition, aircraft_);

        // Step 3: Solve adjoint for each output
        AdjointSolver adjoint(solver_);
        if (!adjoint.IsValid()) {
            return AeroGradients{};
        }

        // Solve all adjoints at once
        const Index n = cached_solution_.circulation.size();
        MatX rhs(n, 6);  // CL, CD, CY, Cl, Cm, Cn
        rhs.col(0) = force_grads.dCL_dGamma;
        rhs.col(1) = force_grads.dCD_dGamma;
        rhs.col(2) = force_grads.dCY_dGamma;
        rhs.col(3) = force_grads.dCl_dGamma;
        rhs.col(4) = force_grads.dCm_dGamma;
        rhs.col(5) = force_grads.dCn_dGamma;

        MatX lambda = adjoint.SolveAdjointMultiple(rhs);

        // Step 4: Compute RHS sensitivity ∂b/∂param
        VecX db_dparam = ComputeRHSSensitivity(condition, param_index);

        // Step 5: Compute implicit terms λᵀ·(∂b/∂param)
        AeroGradients implicit_grads;
        implicit_grads.dCL = lambda.col(0).dot(db_dparam);
        implicit_grads.dCD = lambda.col(1).dot(db_dparam);
        implicit_grads.dCY = lambda.col(2).dot(db_dparam);
        implicit_grads.dCl = lambda.col(3).dot(db_dparam);
        implicit_grads.dCm = lambda.col(4).dot(db_dparam);
        implicit_grads.dCn = lambda.col(5).dot(db_dparam);

        // Step 6: Compute explicit terms ∂J/∂param|_Γ
        AeroGradients explicit_grads = ComputeExplicitDerivatives(param_index);

        // Step 7: Total derivative = explicit + implicit
        AeroGradients total;
        total.dCL = explicit_grads.dCL + implicit_grads.dCL;
        total.dCD = explicit_grads.dCD + implicit_grads.dCD;
        total.dCY = explicit_grads.dCY + implicit_grads.dCY;
        total.dCl = explicit_grads.dCl + implicit_grads.dCl;
        total.dCm = explicit_grads.dCm + implicit_grads.dCm;
        total.dCn = explicit_grads.dCn + implicit_grads.dCn;
        total.dCDi = explicit_grads.dCDi;

        return total;
    }
};

// ============================================================================
// GradientCalculator Public Interface
// ============================================================================

GradientCalculator::GradientCalculator() : impl_(std::make_unique<Impl>()) {}

GradientCalculator::~GradientCalculator() = default;

void GradientCalculator::SetAircraft(const Aircraft& aircraft) {
    impl_->aircraft_ = aircraft;
    impl_->solver_.SetAircraft(aircraft);
    impl_->aircraft_set_ = true;
    impl_->solution_cached_ = false;
}

void GradientCalculator::SetReferencePoint(const Vec3& ref) {
    impl_->solver_.SetReferencePoint(ref);
}

void GradientCalculator::SetFiniteDiffStep(Scalar) {
    // No-op: analytical derivatives don't use finite differences
}

// ============================================================================
// Derivative Computation Methods
// ============================================================================

AeroGradients GradientCalculator::ComputeAlphaDerivatives(
    const FlightCondition& condition) {
    return impl_->AdjointDerivatives(condition, 0);
}

AeroGradients GradientCalculator::ComputeBetaDerivatives(
    const FlightCondition& condition) {
    return impl_->AdjointDerivatives(condition, 1);
}

AeroGradients GradientCalculator::ComputeRollRateDerivatives(
    const FlightCondition& condition) {
    return impl_->AdjointDerivatives(condition, 2);
}

AeroGradients GradientCalculator::ComputePitchRateDerivatives(
    const FlightCondition& condition) {
    return impl_->AdjointDerivatives(condition, 3);
}

AeroGradients GradientCalculator::ComputeYawRateDerivatives(
    const FlightCondition& condition) {
    return impl_->AdjointDerivatives(condition, 4);
}

AeroGradients GradientCalculator::ComputeMachDerivatives(
    const FlightCondition& condition) {
    // Mach changes the Prandtl-Glauert factor, which affects the AIC matrix.
    // This requires re-solving the system at perturbed Mach numbers.
    // TODO: Implement analytical ∂A/∂M via adjoint with PG chain rule
    const Scalar h = 0.001;

    FlightCondition cond_plus = condition;
    FlightCondition cond_minus = condition;
    cond_plus.mach = condition.mach + h;
    cond_minus.mach = condition.mach - h;

    if (cond_plus.mach >= 1.0) {
        cond_plus.mach = 0.999;
    }
    if (cond_minus.mach < 0.0) {
        cond_minus.mach = 0.0;
    }

    AerodynamicResults res_plus = impl_->ComputeResults(cond_plus);
    AerodynamicResults res_minus = impl_->ComputeResults(cond_minus);

    Scalar inv_2h = 0.5 / h;

    AeroGradients grads;
    grads.dCL = (res_plus.CL - res_minus.CL) * inv_2h;
    grads.dCD = (res_plus.CD - res_minus.CD) * inv_2h;
    grads.dCY = (res_plus.CY - res_minus.CY) * inv_2h;
    grads.dCl = (res_plus.Cl - res_minus.Cl) * inv_2h;
    grads.dCm = (res_plus.Cm - res_minus.Cm) * inv_2h;
    grads.dCn = (res_plus.Cn - res_minus.Cn) * inv_2h;
    grads.dCDi = (res_plus.CDi_ff - res_minus.CDi_ff) * inv_2h;

    return grads;
}

FullGradients GradientCalculator::ComputeFlightConditionDerivatives(
    const FlightCondition& condition) {
    FullGradients grads;

    grads.alpha_derivs = ComputeAlphaDerivatives(condition);
    grads.beta_derivs = ComputeBetaDerivatives(condition);
    grads.p_derivs = ComputeRollRateDerivatives(condition);
    grads.q_derivs = ComputePitchRateDerivatives(condition);
    grads.r_derivs = ComputeYawRateDerivatives(condition);
    grads.mach_derivs = ComputeMachDerivatives(condition);

    return grads;
}

// ============================================================================
// Geometry Gradient Computation
// ============================================================================

namespace {

/// Compute the chord-weighted twist derivative weight.
/// Returns ∂(strip_twist)/∂(section_twist_radians) for either the left or right
/// section. This accounts for the atan2 chord-weighted interpolation used by
/// the mesh generator.
inline double ComputeChordWeightedTwistWeight(double fc,
                                              double section_chord_left,
                                              double section_chord_right,
                                              double section_twist_left,
                                              double section_twist_right,
                                              bool is_left_section) {
    double CL = section_chord_left;
    double CR = section_chord_right;
    double tL = section_twist_left;
    double tR = section_twist_right;

    double sinS = (1.0 - fc) * CL * std::sin(tL) + fc * CR * std::sin(tR);
    double cosS = (1.0 - fc) * CL * std::cos(tL) + fc * CR * std::cos(tR);
    double denom = sinS * sinS + cosS * cosS;

    if (denom < 1e-30)
        return is_left_section ? (1.0 - fc) : fc;

    if (is_left_section) {
        double numer =
            (1.0 - fc) * CL * (std::cos(tL) * cosS + std::sin(tL) * sinS);
        return numer / denom;
    } else {
        double numer = fc * CR * (std::cos(tR) * cosS + std::sin(tR) * sinS);
        return numer / denom;
    }
}

/// Get spacing types for a surface given its mesh index.
/// Handles y-duplication mapping from mesh surface to aircraft surface.
void GetSurfaceSpacing(const Aircraft& aircraft, const VortexLattice& mesh,
                       Index surface_idx, SpacingType& spanwise_spacing,
                       SpacingType& chordwise_spacing) {
    const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
    Index original_surface_idx = surface_info.is_y_duplicated
                                     ? surface_info.original_surface_index
                                     : surface_idx;

    Index aircraft_surface_count = 0;
    for (Index i = 0; i < aircraft.GetNumSurfaces(); ++i) {
        if (aircraft_surface_count == original_surface_idx) {
            spanwise_spacing = aircraft.GetSurface(i).GetSpanwiseSpacing();
            chordwise_spacing = aircraft.GetSurface(i).GetChordwiseSpacing();
            return;
        }
        aircraft_surface_count++;
        if (aircraft.GetSurface(i).IsYDuplicated()) {
            aircraft_surface_count++;
        }
    }
    // Fallback to cosine (the default)
    spanwise_spacing = SpacingType::kCosine;
    chordwise_spacing = SpacingType::kCosine;
}

/// Compute the gradient of the RHS with respect to a geometry parameter.
/// The RHS is b_i = -n_i · V_inf, so ∂b_i/∂θ = -(∂n_i/∂θ) · V_inf
VecX ComputeRHSGeometrySensitivity(const VortexLattice& mesh,
                                   const Aircraft& aircraft,
                                   const FlightCondition& condition,
                                   Index surface_idx, Index section_idx,
                                   GeometryVariable variable) {
    const Index n = mesh.GetNumVortices();
    VecX db_dtheta = VecX::Zero(n);

    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    if (surface_idx >= mesh.GetNumSurfaces())
        return db_dtheta;
    const SurfaceMeshInfo& surf_info = mesh.GetSurfaceInfo(surface_idx);

    // Get spacing types for this surface
    SpacingType spanwise_spacing, chordwise_spacing;
    GetSurfaceSpacing(aircraft, mesh, surface_idx, spanwise_spacing,
                      chordwise_spacing);

    if (variable == GeometryVariable::kTwist) {
        SectionStripMapping mapping =
            BuildSectionStripMapping(aircraft, mesh, surface_idx);
        if (mapping.num_sections == 0)
            return db_dtheta;

        Index num_chordwise = surf_info.num_chordwise;

        auto compute_twist_contribution =
            [&](Index first_strip, Index num_strips, bool is_left_section) {
                // Get section chord and twist for chord-weighted interpolation
                const Strip& first_strip_data = mesh.GetStrip(first_strip);
                const Strip& last_strip_data =
                    mesh.GetStrip(first_strip + num_strips - 1);
                double sec_chord_left = first_strip_data.chord_left;
                double sec_chord_right = last_strip_data.chord_right;
                double sec_twist_left = first_strip_data.twist_angle_left;
                double sec_twist_right = last_strip_data.twist_angle_right;

                for (Index istrip = 0; istrip < num_strips; ++istrip) {
                    Index strip_idx = first_strip + istrip;
                    const Strip& strip = mesh.GetStrip(strip_idx);

                    double fc = ComputeStripInterpolationFraction(
                        istrip, num_strips, spanwise_spacing);
                    double weight = ComputeChordWeightedTwistWeight(
                        fc, sec_chord_left, sec_chord_right, sec_twist_left,
                        sec_twist_right, is_left_section);

                    for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
                        Index vor_idx = strip.first_vortex_index + ichord;
                        const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

                        Vec3 bound_dir =
                            (vor.bound_vortex_right - vor.bound_vortex_left)
                                .normalized();
                        Vec3 n_vec = vor.normal_at_control_point;

                        Vec3 dn_dtwist = bound_dir.cross(n_vec);
                        dn_dtwist *= kDegreesToRadians;
                        dn_dtwist *= weight;

                        Scalar db = -(dn_dtwist.x() * v_inf[0] +
                                      dn_dtwist.y() * v_inf[1] +
                                      dn_dtwist.z() * v_inf[2]);

                        db_dtheta(vor_idx) += db;
                    }
                }
            };

        if (section_idx < mapping.num_sections - 1) {
            auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
            if (num_strips > 0) {
                compute_twist_contribution(first_strip, num_strips, true);
            }
        }

        if (section_idx > 0) {
            auto [first_strip, num_strips] =
                mapping.GetStripRange(section_idx - 1);
            if (num_strips > 0) {
                compute_twist_contribution(first_strip, num_strips, false);
            }
        }
    } else if (variable == GeometryVariable::kChord ||
               variable == GeometryVariable::kXle ||
               variable == GeometryVariable::kYle ||
               variable == GeometryVariable::kZle) {
        SectionStripMapping mapping =
            BuildSectionStripMapping(aircraft, mesh, surface_idx);
        if (mapping.num_sections == 0)
            return db_dtheta;

        Index num_chordwise = surf_info.num_chordwise;

        // Pre-compute chordwise fractions
        std::vector<double> xvr_fracs, xcp_fracs;
        BuildDefaultChordwiseFractions(num_chordwise, xvr_fracs, xcp_fracs,
                                       chordwise_spacing);

        auto compute_interval_contribution = [&](Index first_strip,
                                                 Index num_strips,
                                                 bool is_left_section) {
            for (Index istrip = 0; istrip < num_strips; ++istrip) {
                Index strip_idx = first_strip + istrip;
                const Strip& strip = mesh.GetStrip(strip_idx);

                double f1, f2;
                ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2,
                                          spanwise_spacing);

                double w1 = is_left_section ? (1.0 - f1) : f1;
                double w2 = is_left_section ? (1.0 - f2) : f2;

                for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
                    Index vor_idx = strip.first_vortex_index + ichord;
                    const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

                    double xvr = xvr_fracs[ichord];

                    double d_rv1[3] = {0, 0, 0};
                    double d_rv2[3] = {0, 0, 0};

                    switch (variable) {
                        case GeometryVariable::kChord:
                            d_rv1[0] = w1 * xvr;
                            d_rv2[0] = w2 * xvr;
                            break;
                        case GeometryVariable::kXle:
                            d_rv1[0] = w1;
                            d_rv2[0] = w2;
                            break;
                        case GeometryVariable::kYle:
                            d_rv1[1] = w1;
                            d_rv2[1] = w2;
                            break;
                        case GeometryVariable::kZle:
                            d_rv1[2] = w1;
                            d_rv2[2] = w2;
                            break;
                        default:
                            break;
                    }

                    Vec3 bound_vec =
                        vor.bound_vortex_right - vor.bound_vortex_left;
                    double bound_len = bound_vec.norm();
                    if (bound_len < 1e-12)
                        continue;
                    Vec3 bound_dir = bound_vec / bound_len;

                    double d_rv[3] = {d_rv2[0] - d_rv1[0], d_rv2[1] - d_rv1[1],
                                      d_rv2[2] - d_rv1[2]};

                    double bd_dot_drv = bound_dir.x() * d_rv[0] +
                                        bound_dir.y() * d_rv[1] +
                                        bound_dir.z() * d_rv[2];
                    Vec3 d_bound_dir;
                    d_bound_dir.x() = d_rv[0] / bound_len -
                                      bound_dir.x() * bd_dot_drv / bound_len;
                    d_bound_dir.y() = d_rv[1] / bound_len -
                                      bound_dir.y() * bd_dot_drv / bound_len;
                    d_bound_dir.z() = d_rv[2] / bound_len -
                                      bound_dir.z() * bd_dot_drv / bound_len;

                    Vec3 normal = vor.normal_at_control_point;

                    // Use the actual solver chord direction (NOT bound_dir ×
                    // normal)
                    double ang = strip.twist_angle -
                                 std::atan(vor.camber_slope_at_control_point);
                    Vec3 chord_dir;
                    chord_dir.x() = std::cos(ang);
                    chord_dir.y() = -std::sin(ang) * strip.trefftz_normal_y;
                    chord_dir.z() = -std::sin(ang) * strip.trefftz_normal_z;

                    Vec3 d_normal_unnorm = chord_dir.cross(d_bound_dir);

                    Vec3 n_unnorm = chord_dir.cross(bound_dir);
                    double n_unnorm_len = n_unnorm.norm();
                    if (n_unnorm_len < 1e-12)
                        continue;

                    double n_dot_dnu = normal.x() * d_normal_unnorm.x() +
                                       normal.y() * d_normal_unnorm.y() +
                                       normal.z() * d_normal_unnorm.z();
                    Vec3 d_normal;
                    d_normal.x() = d_normal_unnorm.x() / n_unnorm_len -
                                   normal.x() * n_dot_dnu / n_unnorm_len;
                    d_normal.y() = d_normal_unnorm.y() / n_unnorm_len -
                                   normal.y() * n_dot_dnu / n_unnorm_len;
                    d_normal.z() = d_normal_unnorm.z() / n_unnorm_len -
                                   normal.z() * n_dot_dnu / n_unnorm_len;

                    // Add Trefftz normal derivative contribution for kYle/kZle
                    // When LE positions change, trefftz_normal changes, which
                    // rotates chord_dir, which changes normal
                    if (variable == GeometryVariable::kYle ||
                        variable == GeometryVariable::kZle) {
                        double dy_s = strip.leading_edge_right.y() -
                                      strip.leading_edge_left.y();
                        double dz_s = strip.leading_edge_right.z() -
                                      strip.leading_edge_left.z();
                        double yz_mag = std::sqrt(dy_s * dy_s + dz_s * dz_s);
                        if (yz_mag > 1e-12) {
                            double d_dy = (variable == GeometryVariable::kYle)
                                              ? (w2 - w1)
                                              : 0.0;
                            double d_dz = (variable == GeometryVariable::kZle)
                                              ? (w2 - w1)
                                              : 0.0;
                            double cross_term = dz_s * d_dy - dy_s * d_dz;
                            double yz_mag3 = yz_mag * yz_mag * yz_mag;
                            double d_es_y = dy_s * cross_term / yz_mag3;
                            double d_es_z = dz_s * cross_term / yz_mag3;

                            double angle =
                                strip.twist_angle -
                                std::atan(vor.camber_slope_at_control_point);
                            double sin_a = std::sin(angle);

                            if (std::abs(sin_a) > 1e-14 &&
                                (std::abs(d_es_y) > 1e-30 ||
                                 std::abs(d_es_z) > 1e-30)) {
                                double d_ec[3] = {0.0, -sin_a * d_es_y,
                                                  -sin_a * d_es_z};
                                double eb[3] = {bound_dir.x(), bound_dir.y(),
                                                bound_dir.z()};
                                double d_ecxb[3] = {
                                    d_ec[1] * eb[2] - d_ec[2] * eb[1],
                                    d_ec[2] * eb[0] - d_ec[0] * eb[2],
                                    d_ec[0] * eb[1] - d_ec[1] * eb[0]};

                                // Compute emag = |ec × eb|
                                double cos_a = std::cos(angle);
                                double es_y = strip.trefftz_normal_y;
                                double es_z = strip.trefftz_normal_z;
                                double ec_full[3] = {cos_a, -sin_a * es_y,
                                                     -sin_a * es_z};
                                double ecxb_full[3] = {
                                    ec_full[1] * eb[2] - ec_full[2] * eb[1],
                                    ec_full[2] * eb[0] - ec_full[0] * eb[2],
                                    ec_full[0] * eb[1] - ec_full[1] * eb[0]};
                                double emag =
                                    std::sqrt(ecxb_full[0] * ecxb_full[0] +
                                              ecxb_full[1] * ecxb_full[1] +
                                              ecxb_full[2] * ecxb_full[2]);

                                if (emag > 1e-12) {
                                    double n_dot_decxb =
                                        normal.x() * d_ecxb[0] +
                                        normal.y() * d_ecxb[1] +
                                        normal.z() * d_ecxb[2];
                                    d_normal.x() +=
                                        (d_ecxb[0] - normal.x() * n_dot_decxb) /
                                        emag;
                                    d_normal.y() +=
                                        (d_ecxb[1] - normal.y() * n_dot_decxb) /
                                        emag;
                                    d_normal.z() +=
                                        (d_ecxb[2] - normal.z() * n_dot_decxb) /
                                        emag;
                                }
                            }
                        }
                    }

                    Scalar db =
                        -(d_normal.x() * v_inf[0] + d_normal.y() * v_inf[1] +
                          d_normal.z() * v_inf[2]);

                    db_dtheta(vor_idx) += db;
                }
            }
        };

        if (section_idx < mapping.num_sections - 1) {
            auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
            if (num_strips > 0) {
                compute_interval_contribution(first_strip, num_strips, true);
            }
        }

        if (section_idx > 0) {
            auto [first_strip, num_strips] =
                mapping.GetStripRange(section_idx - 1);
            if (num_strips > 0) {
                compute_interval_contribution(first_strip, num_strips, false);
            }
        }
    }

    return db_dtheta;
}

/// Compute λᵀ·(∂A/∂θ)·Γ for a geometry parameter using analytical AIC
/// gradients.
Scalar ComputeAICGeometrySensitivityForInterval(
    const VortexLattice& mesh, const MeshArrays& arrays, const VecX& lambda,
    const VecX& gamma, Index first_strip, Index num_strips, Index num_chordwise,
    double core_chord_frac, double core_width_frac, double beta,
    GeometryVariable variable, bool is_left_section,
    SpacingType spanwise_spacing, SpacingType chordwise_spacing) {
    Scalar result = 0.0;
    const Index n = arrays.num_vortices;

    // Pre-compute chordwise fractions matching the mesh generator
    std::vector<double> xvr_fracs, xcp_fracs;
    BuildDefaultChordwiseFractions(num_chordwise, xvr_fracs, xcp_fracs,
                                   chordwise_spacing);

    // Get section chord and twist for chord-weighted twist interpolation
    const Strip& first_strip_data = mesh.GetStrip(first_strip);
    const Strip& last_strip_data = mesh.GetStrip(first_strip + num_strips - 1);
    double sec_chord_left = first_strip_data.chord_left;
    double sec_chord_right = last_strip_data.chord_right;
    double sec_twist_left = first_strip_data.twist_angle_left;
    double sec_twist_right = last_strip_data.twist_angle_right;

    for (Index istrip = 0; istrip < num_strips; ++istrip) {
        Index strip_idx = first_strip + istrip;
        const Strip& strip = mesh.GetStrip(strip_idx);

        double f1, f2;
        ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2,
                                  spanwise_spacing);
        double fc = ComputeStripInterpolationFraction(istrip, num_strips,
                                                      spanwise_spacing);

        double w1 = is_left_section ? (1.0 - f1) : f1;
        double w2 = is_left_section ? (1.0 - f2) : f2;
        double wc = is_left_section ? (1.0 - fc) : fc;

        for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
            Index vor_idx = strip.first_vortex_index + ichord;

            double xvr = xvr_fracs[ichord];
            double xcp = xcp_fracs[ichord];

            double d_rv1[3] = {0, 0, 0};
            double d_rv2[3] = {0, 0, 0};
            double d_rc[3] = {0, 0, 0};
            double d_normal[3] = {0, 0, 0};

            switch (variable) {
                case GeometryVariable::kChord:
                    d_rv1[0] = w1 * xvr;
                    d_rv2[0] = w2 * xvr;
                    d_rc[0] = wc * xcp;
                    break;

                case GeometryVariable::kXle:
                    d_rv1[0] = w1;
                    d_rv2[0] = w2;
                    d_rc[0] = wc;
                    break;

                case GeometryVariable::kYle:
                    d_rv1[1] = w1;
                    d_rv2[1] = w2;
                    d_rc[1] = wc;
                    break;

                case GeometryVariable::kZle:
                    d_rv1[2] = w1;
                    d_rv2[2] = w2;
                    d_rc[2] = wc;
                    break;

                case GeometryVariable::kTwist: {
                    const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

                    double dx_b =
                        vor.bound_vortex_right.x() - vor.bound_vortex_left.x();
                    double dy_b =
                        vor.bound_vortex_right.y() - vor.bound_vortex_left.y();
                    double dz_b =
                        vor.bound_vortex_right.z() - vor.bound_vortex_left.z();
                    double bound_len =
                        std::sqrt(dx_b * dx_b + dy_b * dy_b + dz_b * dz_b);
                    double bound_dir[3] = {0, 1, 0};
                    if (bound_len > 1e-12) {
                        bound_dir[0] = dx_b / bound_len;
                        bound_dir[1] = dy_b / bound_len;
                        bound_dir[2] = dz_b / bound_len;
                    }

                    double twist = strip.twist_angle;
                    double camber = vor.camber_slope_at_control_point;
                    double angle = twist - std::atan(camber);
                    double ca = std::cos(angle);
                    double sa = std::sin(angle);

                    double dec[3];
                    dec[0] = -sa;
                    dec[1] = -ca * strip.trefftz_normal_y;
                    dec[2] = -ca * strip.trefftz_normal_z;

                    double ec[3] = {ca, -sa * strip.trefftz_normal_y,
                                    -sa * strip.trefftz_normal_z};
                    double ecxb[3];
                    ecxb[0] = ec[1] * bound_dir[2] - ec[2] * bound_dir[1];
                    ecxb[1] = ec[2] * bound_dir[0] - ec[0] * bound_dir[2];
                    ecxb[2] = ec[0] * bound_dir[1] - ec[1] * bound_dir[0];
                    double emag =
                        std::sqrt(ecxb[0] * ecxb[0] + ecxb[1] * ecxb[1] +
                                  ecxb[2] * ecxb[2]);

                    if (emag > 1e-12) {
                        double decxb[3];
                        decxb[0] =
                            dec[1] * bound_dir[2] - dec[2] * bound_dir[1];
                        decxb[1] =
                            dec[2] * bound_dir[0] - dec[0] * bound_dir[2];
                        decxb[2] =
                            dec[0] * bound_dir[1] - dec[1] * bound_dir[0];

                        double normal[3] = {ecxb[0] / emag, ecxb[1] / emag,
                                            ecxb[2] / emag};
                        double demag = normal[0] * decxb[0] +
                                       normal[1] * decxb[1] +
                                       normal[2] * decxb[2];

                        // Use chord-weighted twist derivative with
                        // kDegreesToRadians (section twist parameter is in
                        // degrees)
                        double twist_weight = ComputeChordWeightedTwistWeight(
                            fc, sec_chord_left, sec_chord_right, sec_twist_left,
                            sec_twist_right, is_left_section);
                        double twist_scale = twist_weight * kDegreesToRadians;

                        d_normal[0] =
                            twist_scale * (decxb[0] - normal[0] * demag) / emag;
                        d_normal[1] =
                            twist_scale * (decxb[1] - normal[1] * demag) / emag;
                        d_normal[2] =
                            twist_scale * (decxb[2] - normal[2] * demag) / emag;
                    }
                    break;
                }
            }

            // For position-affecting variables, compute d_normal from bound
            // direction change
            if (variable != GeometryVariable::kTwist) {
                const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);
                Vec3 bound_vec = vor.bound_vortex_right - vor.bound_vortex_left;
                double bound_len = bound_vec.norm();
                if (bound_len > 1e-12) {
                    Vec3 bound_dir = bound_vec / bound_len;

                    double d_rv[3] = {d_rv2[0] - d_rv1[0], d_rv2[1] - d_rv1[1],
                                      d_rv2[2] - d_rv1[2]};
                    double bd_dot_drv = bound_dir.x() * d_rv[0] +
                                        bound_dir.y() * d_rv[1] +
                                        bound_dir.z() * d_rv[2];

                    Vec3 d_bound_dir;
                    d_bound_dir.x() = d_rv[0] / bound_len -
                                      bound_dir.x() * bd_dot_drv / bound_len;
                    d_bound_dir.y() = d_rv[1] / bound_len -
                                      bound_dir.y() * bd_dot_drv / bound_len;
                    d_bound_dir.z() = d_rv[2] / bound_len -
                                      bound_dir.z() * bd_dot_drv / bound_len;

                    Vec3 normal_vec = vor.normal_at_control_point;

                    // Use the actual solver chord direction (NOT bound_dir ×
                    // normal)
                    double ang = strip.twist_angle -
                                 std::atan(vor.camber_slope_at_control_point);
                    Vec3 chord_dir;
                    chord_dir.x() = std::cos(ang);
                    chord_dir.y() = -std::sin(ang) * strip.trefftz_normal_y;
                    chord_dir.z() = -std::sin(ang) * strip.trefftz_normal_z;

                    Vec3 d_normal_unnorm = chord_dir.cross(d_bound_dir);
                    Vec3 n_unnorm = chord_dir.cross(bound_dir);
                    double n_unnorm_len = n_unnorm.norm();
                    if (n_unnorm_len > 1e-12) {
                        double n_dot_dnu =
                            normal_vec.x() * d_normal_unnorm.x() +
                            normal_vec.y() * d_normal_unnorm.y() +
                            normal_vec.z() * d_normal_unnorm.z();
                        d_normal[0] = d_normal_unnorm.x() / n_unnorm_len -
                                      normal_vec.x() * n_dot_dnu / n_unnorm_len;
                        d_normal[1] = d_normal_unnorm.y() / n_unnorm_len -
                                      normal_vec.y() * n_dot_dnu / n_unnorm_len;
                        d_normal[2] = d_normal_unnorm.z() / n_unnorm_len -
                                      normal_vec.z() * n_dot_dnu / n_unnorm_len;
                    }

                    // Add Trefftz normal derivative for kYle/kZle
                    if (variable == GeometryVariable::kYle ||
                        variable == GeometryVariable::kZle) {
                        double dy_s = strip.leading_edge_right.y() -
                                      strip.leading_edge_left.y();
                        double dz_s = strip.leading_edge_right.z() -
                                      strip.leading_edge_left.z();
                        double yz_mag = std::sqrt(dy_s * dy_s + dz_s * dz_s);
                        if (yz_mag > 1e-12) {
                            double d_dy = (variable == GeometryVariable::kYle)
                                              ? (w2 - w1)
                                              : 0.0;
                            double d_dz = (variable == GeometryVariable::kZle)
                                              ? (w2 - w1)
                                              : 0.0;
                            double cross_term = dz_s * d_dy - dy_s * d_dz;
                            double yz_mag3 = yz_mag * yz_mag * yz_mag;
                            double d_es_y = dy_s * cross_term / yz_mag3;
                            double d_es_z = dz_s * cross_term / yz_mag3;

                            double angle =
                                strip.twist_angle -
                                std::atan(vor.camber_slope_at_control_point);
                            double sin_a = std::sin(angle);

                            if (std::abs(sin_a) > 1e-14 &&
                                (std::abs(d_es_y) > 1e-30 ||
                                 std::abs(d_es_z) > 1e-30)) {
                                double d_ec[3] = {0.0, -sin_a * d_es_y,
                                                  -sin_a * d_es_z};
                                double eb[3] = {bound_dir.x(), bound_dir.y(),
                                                bound_dir.z()};
                                double d_ecxb_t[3] = {
                                    d_ec[1] * eb[2] - d_ec[2] * eb[1],
                                    d_ec[2] * eb[0] - d_ec[0] * eb[2],
                                    d_ec[0] * eb[1] - d_ec[1] * eb[0]};
                                double cos_a = std::cos(angle);
                                double es_y_val = strip.trefftz_normal_y;
                                double es_z_val = strip.trefftz_normal_z;
                                double ec_f[3] = {cos_a, -sin_a * es_y_val,
                                                  -sin_a * es_z_val};
                                double ecxb_f[3] = {
                                    ec_f[1] * eb[2] - ec_f[2] * eb[1],
                                    ec_f[2] * eb[0] - ec_f[0] * eb[2],
                                    ec_f[0] * eb[1] - ec_f[1] * eb[0]};
                                double emag_t =
                                    std::sqrt(ecxb_f[0] * ecxb_f[0] +
                                              ecxb_f[1] * ecxb_f[1] +
                                              ecxb_f[2] * ecxb_f[2]);
                                if (emag_t > 1e-12) {
                                    double nv[3] = {normal_vec.x(),
                                                    normal_vec.y(),
                                                    normal_vec.z()};
                                    double n_dot_d = nv[0] * d_ecxb_t[0] +
                                                     nv[1] * d_ecxb_t[1] +
                                                     nv[2] * d_ecxb_t[2];
                                    d_normal[0] +=
                                        (d_ecxb_t[0] - nv[0] * n_dot_d) /
                                        emag_t;
                                    d_normal[1] +=
                                        (d_ecxb_t[1] - nv[1] * n_dot_d) /
                                        emag_t;
                                    d_normal[2] +=
                                        (d_ecxb_t[2] - nv[2] * n_dot_d) /
                                        emag_t;
                                }
                            }
                        }
                    }
                }
            }

            // Row contribution: d(AIC_row · Γ)/d(rc, normal)
            if (d_rc[0] != 0.0 || d_rc[1] != 0.0 || d_rc[2] != 0.0 ||
                d_normal[0] != 0.0 || d_normal[1] != 0.0 ||
                d_normal[2] != 0.0) {
                double d_rc_out[3] = {0, 0, 0};
                double d_normal_out[3] = {0, 0, 0};

                analytical_aic_row_gradient_rc(
                    arrays.rv1.data(), arrays.rv2.data(),
                    &arrays.rc[3 * vor_idx], &arrays.normal[3 * vor_idx],
                    arrays.chord.data(), arrays.component_idx.data(),
                    static_cast<int>(n), static_cast<int>(vor_idx),
                    arrays.component_idx[vor_idx], core_chord_frac,
                    core_width_frac, beta, gamma.data(), d_rc_out,
                    d_normal_out);

                Scalar row_contrib =
                    d_rc_out[0] * d_rc[0] + d_rc_out[1] * d_rc[1] +
                    d_rc_out[2] * d_rc[2] + d_normal_out[0] * d_normal[0] +
                    d_normal_out[1] * d_normal[1] +
                    d_normal_out[2] * d_normal[2];

                result += lambda(vor_idx) * row_contrib;
            }

            // Column contribution: d(λᵀ·A_col)/d(rv1, rv2)
            if (d_rv1[0] != 0.0 || d_rv1[1] != 0.0 || d_rv1[2] != 0.0 ||
                d_rv2[0] != 0.0 || d_rv2[1] != 0.0 || d_rv2[2] != 0.0) {
                double d_rv1_out[3] = {0, 0, 0};
                double d_rv2_out[3] = {0, 0, 0};

                analytical_aic_column_gradient_rv(
                    arrays.rv1.data(), arrays.rv2.data(), arrays.rc.data(),
                    arrays.normal.data(), arrays.chord.data(),
                    arrays.component_idx.data(), static_cast<int>(n),
                    static_cast<int>(vor_idx), core_chord_frac, core_width_frac,
                    beta, lambda.data(), d_rv1_out, d_rv2_out);

                Scalar col_contrib =
                    d_rv1_out[0] * d_rv1[0] + d_rv1_out[1] * d_rv1[1] +
                    d_rv1_out[2] * d_rv1[2] + d_rv2_out[0] * d_rv2[0] +
                    d_rv2_out[1] * d_rv2[1] + d_rv2_out[2] * d_rv2[2];

                result += gamma(vor_idx) * col_contrib;
            }
        }
    }

    return result;
}

/// Compute AIC sensitivity for a section, considering both adjacent intervals.
Scalar ComputeAICGeometrySensitivity(const VortexLattice& mesh,
                                     const Aircraft& aircraft,
                                     const VecX& lambda, const VecX& gamma,
                                     const FlightCondition& condition,
                                     Index surface_idx, Index section_idx,
                                     GeometryVariable variable) {
    const Index n = mesh.GetNumVortices();
    if (n == 0)
        return 0.0;

    SectionStripMapping mapping =
        BuildSectionStripMapping(aircraft, mesh, surface_idx);
    if (mapping.num_sections == 0)
        return 0.0;

    const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
    Index num_chordwise = surface_info.num_chordwise;

    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));

    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    SpacingType spanwise_spacing, chordwise_spacing;
    GetSurfaceSpacing(aircraft, mesh, surface_idx, spanwise_spacing,
                      chordwise_spacing);

    Scalar result = 0.0;

    if (section_idx < mapping.num_sections - 1) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
        if (num_strips > 0) {
            result += ComputeAICGeometrySensitivityForInterval(
                mesh, arrays, lambda, gamma, first_strip, num_strips,
                num_chordwise, core_chord_frac, core_width_frac, beta, variable,
                true, spanwise_spacing, chordwise_spacing);
        }
    }

    if (section_idx > 0) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
        if (num_strips > 0) {
            result += ComputeAICGeometrySensitivityForInterval(
                mesh, arrays, lambda, gamma, first_strip, num_strips,
                num_chordwise, core_chord_frac, core_width_frac, beta, variable,
                false, spanwise_spacing, chordwise_spacing);
        }
    }

    return result;
}

/// Compute velocity kernel geometry sensitivity for one interval.
Scalar ComputeVelocityKernelGeometrySensitivityForInterval(
    const VortexLattice& mesh, const MeshArrays& arrays,
    const std::vector<double>& dJ_dV, const VecX& gamma, Index first_strip,
    Index num_strips, Index num_chordwise, double core_chord_frac,
    double core_width_frac, double beta, GeometryVariable variable,
    bool is_left_section, SpacingType spanwise_spacing,
    SpacingType chordwise_spacing) {
    Scalar result = 0.0;
    const Index n = arrays.num_vortices;

    // Pre-compute chordwise fractions matching the mesh generator
    std::vector<double> xvr_fracs, xcp_fracs;
    BuildDefaultChordwiseFractions(num_chordwise, xvr_fracs, xcp_fracs,
                                   chordwise_spacing);

    for (Index istrip = 0; istrip < num_strips; ++istrip) {
        Index strip_idx = first_strip + istrip;
        const Strip& strip = mesh.GetStrip(strip_idx);

        double f1, f2;
        ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2,
                                  spanwise_spacing);
        double fc = ComputeStripInterpolationFraction(istrip, num_strips,
                                                      spanwise_spacing);

        double w1 = is_left_section ? (1.0 - f1) : f1;
        double w2 = is_left_section ? (1.0 - f2) : f2;
        double wc = is_left_section ? (1.0 - fc) : fc;

        for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
            Index vor_idx = strip.first_vortex_index + ichord;

            double xvr = xvr_fracs[ichord];
            double xcp = xcp_fracs[ichord];

            double d_rv1[3] = {0, 0, 0};
            double d_rv2[3] = {0, 0, 0};
            double d_rc[3] = {0, 0, 0};

            switch (variable) {
                case GeometryVariable::kChord:
                    d_rv1[0] = w1 * xvr;
                    d_rv2[0] = w2 * xvr;
                    d_rc[0] = wc * xcp;
                    break;
                case GeometryVariable::kXle:
                    d_rv1[0] = w1;
                    d_rv2[0] = w2;
                    d_rc[0] = wc;
                    break;
                case GeometryVariable::kYle:
                    d_rv1[1] = w1;
                    d_rv2[1] = w2;
                    d_rc[1] = wc;
                    break;
                case GeometryVariable::kZle:
                    d_rv1[2] = w1;
                    d_rv2[2] = w2;
                    d_rc[2] = wc;
                    break;
                case GeometryVariable::kTwist:
                    break;
            }

            // Row contribution: when control point moves
            if (d_rc[0] != 0.0 || d_rc[1] != 0.0 || d_rc[2] != 0.0) {
                double d_rc_out[3] = {0, 0, 0};
                const double* dJ_dV_i = &dJ_dV[3 * vor_idx];

                analytical_velocity_row_gradient_rc(
                    arrays.rv1.data(), arrays.rv2.data(),
                    &arrays.rc[3 * vor_idx], arrays.chord.data(),
                    arrays.component_idx.data(), static_cast<int>(n),
                    static_cast<int>(vor_idx), arrays.component_idx[vor_idx],
                    core_chord_frac, core_width_frac, beta, gamma.data(),
                    dJ_dV_i, d_rc_out);

                result += d_rc_out[0] * d_rc[0] + d_rc_out[1] * d_rc[1] +
                          d_rc_out[2] * d_rc[2];
            }

            // Column contribution: when bound vortex moves
            if (d_rv1[0] != 0.0 || d_rv1[1] != 0.0 || d_rv1[2] != 0.0 ||
                d_rv2[0] != 0.0 || d_rv2[1] != 0.0 || d_rv2[2] != 0.0) {
                double d_rv1_out[3] = {0, 0, 0};
                double d_rv2_out[3] = {0, 0, 0};

                analytical_velocity_column_gradient_rv(
                    arrays.rv1.data(), arrays.rv2.data(), arrays.rc.data(),
                    arrays.chord.data(), arrays.component_idx.data(),
                    static_cast<int>(n), static_cast<int>(vor_idx),
                    core_chord_frac, core_width_frac, beta, gamma(vor_idx),
                    dJ_dV.data(), d_rv1_out, d_rv2_out);

                result += d_rv1_out[0] * d_rv1[0] + d_rv1_out[1] * d_rv1[1] +
                          d_rv1_out[2] * d_rv1[2];
                result += d_rv2_out[0] * d_rv2[0] + d_rv2_out[1] * d_rv2[1] +
                          d_rv2_out[2] * d_rv2[2];
            }
        }
    }

    return result;
}

/// Compute velocity kernel geometry sensitivity for a section.
Scalar ComputeVelocityKernelGeometrySensitivity(
    const VortexLattice& mesh, const Aircraft& aircraft,
    const std::vector<double>& dJ_dV, const VecX& gamma,
    const FlightCondition& condition, Index surface_idx, Index section_idx,
    GeometryVariable variable) {
    SectionStripMapping mapping =
        BuildSectionStripMapping(aircraft, mesh, surface_idx);
    if (mapping.num_sections == 0)
        return 0.0;

    const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
    Index num_chordwise = surface_info.num_chordwise;

    MeshArrays arrays = ExtractMeshArrays(mesh);
    double beta =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    SpacingType spanwise_spacing, chordwise_spacing;
    GetSurfaceSpacing(aircraft, mesh, surface_idx, spanwise_spacing,
                      chordwise_spacing);

    Scalar result = 0.0;

    if (section_idx < mapping.num_sections - 1) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
        if (num_strips > 0) {
            result += ComputeVelocityKernelGeometrySensitivityForInterval(
                mesh, arrays, dJ_dV, gamma, first_strip, num_strips,
                num_chordwise, core_chord_frac, core_width_frac, beta, variable,
                true, spanwise_spacing, chordwise_spacing);
        }
    }

    if (section_idx > 0) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
        if (num_strips > 0) {
            result += ComputeVelocityKernelGeometrySensitivityForInterval(
                mesh, arrays, dJ_dV, gamma, first_strip, num_strips,
                num_chordwise, core_chord_frac, core_width_frac, beta, variable,
                false, spanwise_spacing, chordwise_spacing);
        }
    }

    return result;
}

/// Compute explicit force term ∂J/∂θ|_Γ,V for a geometry parameter.
AeroGradients ComputeExplicitForceTerm(const VortexLattice& mesh,
                                       const Aircraft& aircraft,
                                       const SolutionState& solution,
                                       const FlightCondition& condition,
                                       Index surface_idx, Index section_idx,
                                       GeometryVariable variable) {
    const Index n = mesh.GetNumVortices();
    if (n == 0)
        return AeroGradients{};

    SectionStripMapping mapping =
        BuildSectionStripMapping(aircraft, mesh, surface_idx);
    if (mapping.num_sections == 0)
        return AeroGradients{};

    const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
    Index num_chordwise = surface_info.num_chordwise;

    Vec3 moment_ref = aircraft.GetMomentReferencePoint();
    std::vector<double> bound_left(3 * n);
    std::vector<double> bound_right(3 * n);
    std::vector<double> moment_arm(3 * n);
    std::vector<double> induced_velocity(3 * n);

    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);
        bound_left[3 * i + 0] = vor.bound_vortex_left.x();
        bound_left[3 * i + 1] = vor.bound_vortex_left.y();
        bound_left[3 * i + 2] = vor.bound_vortex_left.z();
        bound_right[3 * i + 0] = vor.bound_vortex_right.x();
        bound_right[3 * i + 1] = vor.bound_vortex_right.y();
        bound_right[3 * i + 2] = vor.bound_vortex_right.z();
        moment_arm[3 * i + 0] = vor.control_point.x() - moment_ref.x();
        moment_arm[3 * i + 1] = vor.control_point.y() - moment_ref.y();
        moment_arm[3 * i + 2] = vor.control_point.z() - moment_ref.z();
        if (solution.induced_velocity.rows() > i) {
            induced_velocity[3 * i + 0] = solution.induced_velocity(i, 0);
            induced_velocity[3 * i + 1] = solution.induced_velocity(i, 1);
            induced_velocity[3 * i + 2] = solution.induced_velocity(i, 2);
        } else {
            induced_velocity[3 * i + 0] = 0.0;
            induced_velocity[3 * i + 1] = 0.0;
            induced_velocity[3 * i + 2] = 0.0;
        }
    }

    std::vector<double> circulation(n);
    for (Index i = 0; i < n; ++i) {
        circulation[i] = solution.circulation(i);
    }

    // Get spacing types
    SpacingType spanwise_spacing, chordwise_spacing;
    GetSurfaceSpacing(aircraft, mesh, surface_idx, spanwise_spacing,
                      chordwise_spacing);

    // Pre-compute chordwise fractions
    std::vector<double> xvr_fracs, xcp_fracs;
    BuildDefaultChordwiseFractions(num_chordwise, xvr_fracs, xcp_fracs,
                                   chordwise_spacing);

    // Initialize perturbation arrays to zero
    std::vector<double> d_bound_left(3 * n, 0.0);
    std::vector<double> d_bound_right(3 * n, 0.0);
    std::vector<double> d_moment_arm(3 * n, 0.0);

    auto compute_interval_perturbations =
        [&](Index first_strip, Index num_strips, bool is_left_section) {
            for (Index istrip = 0; istrip < num_strips; ++istrip) {
                Index strip_idx = first_strip + istrip;
                const Strip& strip = mesh.GetStrip(strip_idx);

                double f1, f2;
                ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2,
                                          spanwise_spacing);
                double fc = ComputeStripInterpolationFraction(
                    istrip, num_strips, spanwise_spacing);

                double w1 = is_left_section ? (1.0 - f1) : f1;
                double w2 = is_left_section ? (1.0 - f2) : f2;
                double wc = is_left_section ? (1.0 - fc) : fc;

                for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
                    Index vor_idx = strip.first_vortex_index + ichord;

                    double xvr = xvr_fracs[ichord];
                    double xcp = xcp_fracs[ichord];

                    switch (variable) {
                        case GeometryVariable::kChord:
                            d_bound_left[3 * vor_idx + 0] += w1 * xvr;
                            d_bound_right[3 * vor_idx + 0] += w2 * xvr;
                            d_moment_arm[3 * vor_idx + 0] += wc * xcp;
                            break;

                        case GeometryVariable::kXle:
                            d_bound_left[3 * vor_idx + 0] += w1;
                            d_bound_right[3 * vor_idx + 0] += w2;
                            d_moment_arm[3 * vor_idx + 0] += wc;
                            break;

                        case GeometryVariable::kYle:
                            d_bound_left[3 * vor_idx + 1] += w1;
                            d_bound_right[3 * vor_idx + 1] += w2;
                            d_moment_arm[3 * vor_idx + 1] += wc;
                            break;

                        case GeometryVariable::kZle:
                            d_bound_left[3 * vor_idx + 2] += w1;
                            d_bound_right[3 * vor_idx + 2] += w2;
                            d_moment_arm[3 * vor_idx + 2] += wc;
                            break;

                        case GeometryVariable::kTwist:
                            break;
                    }
                }
            }
        };

    if (section_idx < mapping.num_sections - 1) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);
        if (num_strips > 0) {
            compute_interval_perturbations(first_strip, num_strips, true);
        }
    }

    if (section_idx > 0) {
        auto [first_strip, num_strips] = mapping.GetStripRange(section_idx - 1);
        if (num_strips > 0) {
            compute_interval_perturbations(first_strip, num_strips, false);
        }
    }

    double d_result[6] = {0};
    double ref_area = aircraft.GetReferenceArea();

    analytical_force_geometry_gradient(
        condition.alpha, condition.beta, circulation.data(),
        induced_velocity.data(), bound_left.data(), bound_right.data(),
        moment_arm.data(), static_cast<int>(n), ref_area, d_bound_left.data(),
        d_bound_right.data(), d_moment_arm.data(), d_result);

    AeroGradients grads;
    grads.dCL = d_result[0];
    grads.dCD = d_result[1];
    grads.dCY = d_result[2];
    grads.dCl = d_result[3];
    grads.dCm = d_result[4];
    grads.dCn = d_result[5];
    grads.dCDi = 0.0;

    return grads;
}

/// Compute geometry gradient using the adjoint method with analytical
/// derivatives.
AeroGradients ComputeSectionGradientAdjoint(const Aircraft& aircraft,
                                            const FlightCondition& condition,
                                            Index surface_idx,
                                            Index section_idx,
                                            GeometryVariable variable) {
    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);
    SolutionState solution = solver.Solve(condition);

    const VortexLattice& mesh = solver.GetMesh();
    const Index n = solution.circulation.size();

    if (n == 0) {
        return AeroGradients{};
    }

    // Step 1: Compute force gradients ∂J/∂Γ
    ForceGradients force_grads =
        ComputeForceGradients(mesh, solution, condition, aircraft);

    // Step 2: Solve adjoint system Aᵀλ = ∂J/∂Γ
    AdjointSolver adjoint(solver);
    if (!adjoint.IsValid()) {
        return AeroGradients{};
    }

    MatX rhs(n, 6);
    rhs.col(0) = force_grads.dCL_dGamma;
    rhs.col(1) = force_grads.dCD_dGamma;
    rhs.col(2) = force_grads.dCY_dGamma;
    rhs.col(3) = force_grads.dCl_dGamma;
    rhs.col(4) = force_grads.dCm_dGamma;
    rhs.col(5) = force_grads.dCn_dGamma;

    MatX lambda = adjoint.SolveAdjointMultiple(rhs);

    // Step 3: Compute RHS sensitivity ∂b/∂θ
    VecX db_dtheta = ComputeRHSGeometrySensitivity(
        mesh, aircraft, condition, surface_idx, section_idx, variable);

    // Step 4: Compute implicit term λᵀ·(∂b/∂θ)
    AeroGradients implicit_rhs;
    implicit_rhs.dCL = lambda.col(0).dot(db_dtheta);
    implicit_rhs.dCD = lambda.col(1).dot(db_dtheta);
    implicit_rhs.dCY = lambda.col(2).dot(db_dtheta);
    implicit_rhs.dCl = lambda.col(3).dot(db_dtheta);
    implicit_rhs.dCm = lambda.col(4).dot(db_dtheta);
    implicit_rhs.dCn = lambda.col(5).dot(db_dtheta);
    implicit_rhs.dCDi = 0.0;

    // Step 5: Compute AIC sensitivity term -λᵀ·(∂A/∂θ)·Γ
    AeroGradients aic_sens;
    aic_sens.dCL = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(0), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCD = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(1), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCY = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(2), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCl = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(3), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCm = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(4), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCn = -ComputeAICGeometrySensitivity(
        mesh, aircraft, lambda.col(5), solution.circulation, condition,
        surface_idx, section_idx, variable);
    aic_sens.dCDi = 0.0;

    // Step 6: Compute explicit term ∂J/∂θ|_Γ,V
    AeroGradients explicit_term =
        ComputeExplicitForceTerm(mesh, aircraft, solution, condition,
                                 surface_idx, section_idx, variable);

    // Step 7: Compute velocity kernel geometry sensitivity (∂J/∂V)·(∂K/∂θ)·Γ
    AeroGradients velocity_kernel_sens;
    velocity_kernel_sens.dCL = 0.0;
    velocity_kernel_sens.dCD = 0.0;
    velocity_kernel_sens.dCY = 0.0;
    velocity_kernel_sens.dCl = 0.0;
    velocity_kernel_sens.dCm = 0.0;
    velocity_kernel_sens.dCn = 0.0;
    velocity_kernel_sens.dCDi = 0.0;

    {
        Vec3 moment_ref = aircraft.GetMomentReferencePoint();
        std::vector<double> bound_left(3 * n);
        std::vector<double> bound_right(3 * n);
        std::vector<double> moment_arm(3 * n);
        std::vector<double> induced_velocity(3 * n);
        std::vector<double> circulation(n);

        for (Index i = 0; i < n; ++i) {
            const HorseshoeVortex& vor = mesh.GetVortex(i);
            bound_left[3 * i + 0] = vor.bound_vortex_left.x();
            bound_left[3 * i + 1] = vor.bound_vortex_left.y();
            bound_left[3 * i + 2] = vor.bound_vortex_left.z();
            bound_right[3 * i + 0] = vor.bound_vortex_right.x();
            bound_right[3 * i + 1] = vor.bound_vortex_right.y();
            bound_right[3 * i + 2] = vor.bound_vortex_right.z();
            moment_arm[3 * i + 0] = vor.control_point.x() - moment_ref.x();
            moment_arm[3 * i + 1] = vor.control_point.y() - moment_ref.y();
            moment_arm[3 * i + 2] = vor.control_point.z() - moment_ref.z();
            if (solution.induced_velocity.rows() > i) {
                induced_velocity[3 * i + 0] = solution.induced_velocity(i, 0);
                induced_velocity[3 * i + 1] = solution.induced_velocity(i, 1);
                induced_velocity[3 * i + 2] = solution.induced_velocity(i, 2);
            }
            circulation[i] = solution.circulation(i);
        }

        double ref_area = aircraft.GetReferenceArea();

        auto compute_velocity_kernel_term = [&](int output_index) -> Scalar {
            std::vector<double> dJ_dV(3 * n, 0.0);
            analytical_force_induced_velocity_gradient(
                condition.alpha, condition.beta, circulation.data(),
                induced_velocity.data(), bound_left.data(), bound_right.data(),
                nullptr, moment_arm.data(), static_cast<int>(n), ref_area,
                output_index, dJ_dV.data());

            return ComputeVelocityKernelGeometrySensitivity(
                mesh, aircraft, dJ_dV, solution.circulation, condition,
                surface_idx, section_idx, variable);
        };

        velocity_kernel_sens.dCL = compute_velocity_kernel_term(0);
        velocity_kernel_sens.dCD = compute_velocity_kernel_term(1);
        velocity_kernel_sens.dCY = compute_velocity_kernel_term(2);
        velocity_kernel_sens.dCl = compute_velocity_kernel_term(3);
        velocity_kernel_sens.dCm = compute_velocity_kernel_term(4);
        velocity_kernel_sens.dCn = compute_velocity_kernel_term(5);
    }

    // Total gradient = explicit + velocity_kernel + implicit_rhs + aic_sens
    AeroGradients result;
    result.dCL = explicit_term.dCL + velocity_kernel_sens.dCL +
                 implicit_rhs.dCL + aic_sens.dCL;
    result.dCD = explicit_term.dCD + velocity_kernel_sens.dCD +
                 implicit_rhs.dCD + aic_sens.dCD;
    result.dCY = explicit_term.dCY + velocity_kernel_sens.dCY +
                 implicit_rhs.dCY + aic_sens.dCY;
    result.dCl = explicit_term.dCl + velocity_kernel_sens.dCl +
                 implicit_rhs.dCl + aic_sens.dCl;
    result.dCm = explicit_term.dCm + velocity_kernel_sens.dCm +
                 implicit_rhs.dCm + aic_sens.dCm;
    result.dCn = explicit_term.dCn + velocity_kernel_sens.dCn +
                 implicit_rhs.dCn + aic_sens.dCn;
    result.dCDi = 0.0;  // TODO: Add Trefftz plane gradient

#ifndef NDEBUG
    if (std::getenv("VLM_DEBUG_ADJOINT")) {
        fprintf(stderr, "Adjoint terms for var=%d sec=%lld:\n",
                static_cast<int>(variable),
                static_cast<long long>(section_idx));
        fprintf(stderr, "  explicit:      dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
                explicit_term.dCL, explicit_term.dCD, explicit_term.dCm);
        fprintf(stderr, "  vel_kernel:    dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
                velocity_kernel_sens.dCL, velocity_kernel_sens.dCD,
                velocity_kernel_sens.dCm);
        fprintf(stderr, "  implicit_rhs:  dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
                implicit_rhs.dCL, implicit_rhs.dCD, implicit_rhs.dCm);
        fprintf(stderr, "  aic_sens:      dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
                aic_sens.dCL, aic_sens.dCD, aic_sens.dCm);
        fprintf(stderr, "  TOTAL:         dCL=%+.6f dCD=%+.6f dCm=%+.6f\n",
                result.dCL, result.dCD, result.dCm);
    }
#endif

    return result;
}

}  // namespace

AeroGradients GradientCalculator::ComputeSectionGradient(
    const FlightCondition& condition, Index surface_index, Index section_index,
    GeometryVariable variable) {
    if (!impl_->aircraft_set_) {
        return AeroGradients{};
    }

    return ComputeSectionGradientAdjoint(
        impl_->aircraft_, condition, surface_index, section_index, variable);
}

GeometryGradients GradientCalculator::ComputeGeometryGradients(
    const FlightCondition& condition) {
    GeometryGradients grads;

    if (!impl_->aircraft_set_) {
        return grads;
    }

    const Aircraft& aircraft = impl_->aircraft_;

    for (Index surf_idx = 0;
         surf_idx < static_cast<Index>(aircraft.GetNumSurfaces()); ++surf_idx) {
        const Surface& surface = aircraft.GetSurface(surf_idx);

        for (Index sect_idx = 0;
             sect_idx < static_cast<Index>(surface.GetNumSections());
             ++sect_idx) {
            SectionGeometryGradients section_grads;
            section_grads.surface_name = surface.GetName();
            section_grads.section_index = sect_idx;

            section_grads.chord_derivs = ComputeSectionGradient(
                condition, surf_idx, sect_idx, GeometryVariable::kChord);
            section_grads.twist_derivs = ComputeSectionGradient(
                condition, surf_idx, sect_idx, GeometryVariable::kTwist);
            section_grads.xle_derivs = ComputeSectionGradient(
                condition, surf_idx, sect_idx, GeometryVariable::kXle);
            section_grads.yle_derivs = ComputeSectionGradient(
                condition, surf_idx, sect_idx, GeometryVariable::kYle);
            section_grads.zle_derivs = ComputeSectionGradient(
                condition, surf_idx, sect_idx, GeometryVariable::kZle);

            grads.section_gradients.push_back(section_grads);
        }
    }

    return grads;
}

// ============================================================================
// Mesh Input Methods
// ============================================================================

void GradientCalculator::SetMesh(const MeshInput& input) {
    impl_->mesh_input_ = input;
    impl_->ref_quantities_ = input.reference;
    impl_->solver_.SetMesh(input);
    impl_->mesh_set_ = true;
    impl_->aircraft_set_ = false;
    impl_->solution_cached_ = false;

    // Sync aircraft reference quantities so AdjointDerivatives (which uses
    // aircraft_ for ComputeForceGradients and ComputeExplicitDerivatives)
    // picks up the correct reference area and moment reference point.
    impl_->aircraft_.SetReferenceArea(input.reference.reference_area);
    impl_->aircraft_.SetReferenceChord(input.reference.reference_chord);
    impl_->aircraft_.SetReferenceSpan(input.reference.reference_span);
    impl_->aircraft_.SetMomentReferencePoint(
        input.reference.moment_reference_point);
}

std::vector<SurfaceMeshNodeGradients>
GradientCalculator::ComputeMeshNodeGradients(const FlightCondition& condition) {
    std::vector<SurfaceMeshNodeGradients> result;
    if (!impl_->mesh_set_)
        return result;

    const MeshInput& mesh_input = impl_->mesh_input_;
    const ReferenceQuantities& ref = impl_->ref_quantities_;

    // Step 1: Solve forward problem
    if (!impl_->solution_cached_ ||
        condition.alpha != impl_->cached_condition_.alpha ||
        condition.beta != impl_->cached_condition_.beta) {
        impl_->SolveAndCache(condition);
    }

    const VortexLattice& mesh = impl_->solver_.GetMesh();
    const SolutionState& solution = impl_->cached_solution_;
    const Index n = solution.circulation.size();
    if (n == 0)
        return result;

    // Step 2: Compute force gradients dJ/dGamma for all 6 outputs
    ForceGradients force_grads =
        ComputeForceGradients(mesh, solution, condition, ref);

    // Step 3: Solve adjoint systems A^T * lambda = dJ/dGamma
    AdjointSolver adjoint(impl_->solver_);
    if (!adjoint.IsValid())
        return result;

    MatX rhs(n, 6);
    rhs.col(0) = force_grads.dCL_dGamma;
    rhs.col(1) = force_grads.dCD_dGamma;
    rhs.col(2) = force_grads.dCY_dGamma;
    rhs.col(3) = force_grads.dCl_dGamma;
    rhs.col(4) = force_grads.dCm_dGamma;
    rhs.col(5) = force_grads.dCn_dGamma;

    MatX lambda = adjoint.SolveAdjointMultiple(rhs);

    // Extract mesh arrays for analytical kernels
    MeshArrays arrays = ExtractMeshArrays(mesh);

    double beta_pg =
        std::sqrt(std::max(0.01, 1.0 - condition.mach * condition.mach));
    double core_chord_frac = mesh.GetVortexCoreRadiusChordFraction();
    double core_width_frac = mesh.GetVortexCoreRadiusWidthFraction();

    // Freestream direction
    double v_inf[3];
    compute_freestream_direction(condition.alpha, condition.beta, v_inf);

    // Extract geometric data for force kernel
    Vec3 moment_ref = ref.moment_reference_point;
    std::vector<double> bound_left(3 * n), bound_right(3 * n);
    std::vector<double> moment_arm(3 * n), induced_velocity(3 * n);
    std::vector<double> circulation(n);
    for (Index i = 0; i < n; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);
        bound_left[3 * i + 0] = vor.bound_vortex_left.x();
        bound_left[3 * i + 1] = vor.bound_vortex_left.y();
        bound_left[3 * i + 2] = vor.bound_vortex_left.z();
        bound_right[3 * i + 0] = vor.bound_vortex_right.x();
        bound_right[3 * i + 1] = vor.bound_vortex_right.y();
        bound_right[3 * i + 2] = vor.bound_vortex_right.z();
        moment_arm[3 * i + 0] = vor.control_point.x() - moment_ref.x();
        moment_arm[3 * i + 1] = vor.control_point.y() - moment_ref.y();
        moment_arm[3 * i + 2] = vor.control_point.z() - moment_ref.z();
        if (solution.induced_velocity.rows() > i) {
            induced_velocity[3 * i + 0] = solution.induced_velocity(i, 0);
            induced_velocity[3 * i + 1] = solution.induced_velocity(i, 1);
            induced_velocity[3 * i + 2] = solution.induced_velocity(i, 2);
        }
        circulation[i] = solution.circulation(i);
    }

    // Step 4: For each objective, accumulate full MeshGradients then chain-rule
    // to nodes We compute for all 6 objectives
    for (int obj = 0; obj < 6; ++obj) {
        VecX lam = lambda.col(obj);

        MeshGradients mesh_grads;
        mesh_grads.Resize(n);

        // (a) RHS term: lambda^T * (db/d_mesh)
        // b_i = -n_i . V_inf, so db/d(normal_i) = -V_inf
        // Contribution: d_normal[i] += -lambda[i] * V_inf (reversed sign in
        // chain rule) But the total derivative has +lambda^T * db/dtheta, so:
        // d_normal[i] += lambda[i] * (-V_inf)  ... wait, let me think.
        //
        // Total: dJ/dtheta = explicit + lambda^T * db/dtheta - lambda^T *
        // dA/dtheta * Gamma We accumulate into MeshGradients the "seed" that
        // when contracted with d(mesh_positions)/d(mesh_nodes), gives the total
        // derivative.
        //
        // For the RHS term: lambda^T * db/dtheta
        //   b_i = -normal_i . V_inf (only normal depends on mesh geometry)
        //   db_i/d(normal_i) = -V_inf
        //   So lambda^T * db/d(normal) contributes: d_normal[i] += lambda[i] *
        //   (-V_inf[k]) for each k
        for (Index i = 0; i < n; ++i) {
            mesh_grads.d_normal[3 * i + 0] += lam(i) * (-v_inf[0]);
            mesh_grads.d_normal[3 * i + 1] += lam(i) * (-v_inf[1]);
            mesh_grads.d_normal[3 * i + 2] += lam(i) * (-v_inf[2]);
        }

        // (b) AIC term: -lambda^T * dA/dtheta * Gamma
        // For row i: the AIC row depends on rc_i and normal_i
        // For column j: the AIC column depends on rv1_j and rv2_j
        //
        // Row gradient: d(AIC_row[i] . Gamma)/d(rc_i) and d/d(normal_i)
        // Contribution scaled by -lambda[i]
        for (Index i = 0; i < n; ++i) {
            double d_rc_out[3] = {0, 0, 0};
            double d_normal_out[3] = {0, 0, 0};

            analytical_aic_row_gradient_rc(
                arrays.rv1.data(), arrays.rv2.data(), &arrays.rc[3 * i],
                &arrays.normal[3 * i], arrays.chord.data(),
                arrays.component_idx.data(), static_cast<int>(n),
                static_cast<int>(i), arrays.component_idx[i], core_chord_frac,
                core_width_frac, beta_pg, solution.circulation.data(), d_rc_out,
                d_normal_out);

            for (int k = 0; k < 3; ++k) {
                mesh_grads.d_rc[3 * i + k] += -lam(i) * d_rc_out[k];
                mesh_grads.d_normal[3 * i + k] += -lam(i) * d_normal_out[k];
            }
        }

        // Column gradient: d(lambda^T . AIC_col[j])/d(rv1_j, rv2_j)
        // Contribution scaled by -Gamma[j]
        for (Index j = 0; j < n; ++j) {
            double d_rv1_out[3] = {0, 0, 0};
            double d_rv2_out[3] = {0, 0, 0};

            analytical_aic_column_gradient_rv(
                arrays.rv1.data(), arrays.rv2.data(), arrays.rc.data(),
                arrays.normal.data(), arrays.chord.data(),
                arrays.component_idx.data(), static_cast<int>(n),
                static_cast<int>(j), core_chord_frac, core_width_frac, beta_pg,
                lam.data(), d_rv1_out, d_rv2_out);

            for (int k = 0; k < 3; ++k) {
                mesh_grads.d_rv1[3 * j + k] +=
                    -solution.circulation(j) * d_rv1_out[k];
                mesh_grads.d_rv2[3 * j + k] +=
                    -solution.circulation(j) * d_rv2_out[k];
            }
        }

        // (c) Explicit force term: dJ/d(rv1, rv2, moment_arm)|_Gamma,V
        // We use the analytical_force_geometry_gradient with unit perturbations
        // per vortex. But that function is a forward-mode computation that
        // takes d_bound_left, d_bound_right, d_moment_arm perturbation vectors.
        // For the adjoint approach, we need to accumulate contributions per
        // vortex.
        //
        // The force is: F_total = sum_i 2*Gamma_i * V_eff_i x (rv2_i - rv1_i)
        // And moment contributions involve moment_arm.
        // We can compute the explicit per-vortex sensitivities.
        {
            // Per-vortex perturbations: set one at a time for each vortex
            // Actually, this is expensive (O(N^2)). Instead, let's use the
            // vectorized force gradient: pass unit perturbation for each
            // vortex. The force geometry gradient is linear in the
            // perturbation, so we can pass identity perturbation per vortex.
            for (Index vi = 0; vi < n; ++vi) {
                // Perturb rv1 (bound_left) of vortex vi in each direction
                for (int dir = 0; dir < 3; ++dir) {
                    std::vector<double> d_bl(3 * n, 0.0), d_br(3 * n, 0.0),
                        d_ma(3 * n, 0.0);
                    d_bl[3 * vi + dir] = 1.0;
                    double d_res[6];
                    analytical_force_geometry_gradient(
                        condition.alpha, condition.beta, circulation.data(),
                        induced_velocity.data(), bound_left.data(),
                        bound_right.data(), moment_arm.data(),
                        static_cast<int>(n), ref.reference_area, d_bl.data(),
                        d_br.data(), d_ma.data(), d_res);
                    mesh_grads.d_rv1[3 * vi + dir] += d_res[obj];
                }
                // Perturb rv2 (bound_right) of vortex vi in each direction
                for (int dir = 0; dir < 3; ++dir) {
                    std::vector<double> d_bl(3 * n, 0.0), d_br(3 * n, 0.0),
                        d_ma(3 * n, 0.0);
                    d_br[3 * vi + dir] = 1.0;
                    double d_res[6];
                    analytical_force_geometry_gradient(
                        condition.alpha, condition.beta, circulation.data(),
                        induced_velocity.data(), bound_left.data(),
                        bound_right.data(), moment_arm.data(),
                        static_cast<int>(n), ref.reference_area, d_bl.data(),
                        d_br.data(), d_ma.data(), d_res);
                    mesh_grads.d_rv2[3 * vi + dir] += d_res[obj];
                }
                // Perturb moment_arm (maps to rc) of vortex vi in each
                // direction
                for (int dir = 0; dir < 3; ++dir) {
                    std::vector<double> d_bl(3 * n, 0.0), d_br(3 * n, 0.0),
                        d_ma(3 * n, 0.0);
                    d_ma[3 * vi + dir] = 1.0;
                    double d_res[6];
                    analytical_force_geometry_gradient(
                        condition.alpha, condition.beta, circulation.data(),
                        induced_velocity.data(), bound_left.data(),
                        bound_right.data(), moment_arm.data(),
                        static_cast<int>(n), ref.reference_area, d_bl.data(),
                        d_br.data(), d_ma.data(), d_res);
                    mesh_grads.d_rc[3 * vi + dir] += d_res[obj];
                }
            }
        }

        // (d) Velocity kernel geometry term: (dJ/dV) . (dV/d_mesh) .
        // d(mesh)/d(nodes) dJ/dV is computed analytically, then we propagate
        // through the velocity kernel
        {
            std::vector<double> dJ_dV(3 * n, 0.0);
            analytical_force_induced_velocity_gradient(
                condition.alpha, condition.beta, circulation.data(),
                induced_velocity.data(), bound_left.data(), bound_right.data(),
                nullptr, moment_arm.data(), static_cast<int>(n),
                ref.reference_area, obj, dJ_dV.data());

            // Propagate through velocity kernel: for each vortex, compute
            // d(sum_i dJ_dV . V[i,j]*Gamma[j])/d(rc_i, rv1_j, rv2_j)
            for (Index i = 0; i < n; ++i) {
                double d_rc_out[3] = {0, 0, 0};
                analytical_velocity_row_gradient_rc(
                    arrays.rv1.data(), arrays.rv2.data(), &arrays.rc[3 * i],
                    arrays.chord.data(), arrays.component_idx.data(),
                    static_cast<int>(n), static_cast<int>(i),
                    arrays.component_idx[i], core_chord_frac, core_width_frac,
                    beta_pg, solution.circulation.data(), &dJ_dV[3 * i],
                    d_rc_out);
                for (int k = 0; k < 3; ++k) {
                    mesh_grads.d_rc[3 * i + k] += d_rc_out[k];
                }
            }

            for (Index j = 0; j < n; ++j) {
                double d_rv1_out[3] = {0, 0, 0};
                double d_rv2_out[3] = {0, 0, 0};
                analytical_velocity_column_gradient_rv(
                    arrays.rv1.data(), arrays.rv2.data(), arrays.rc.data(),
                    arrays.chord.data(), arrays.component_idx.data(),
                    static_cast<int>(n), static_cast<int>(j), core_chord_frac,
                    core_width_frac, beta_pg, solution.circulation(j),
                    dJ_dV.data(), d_rv1_out, d_rv2_out);
                for (int k = 0; k < 3; ++k) {
                    mesh_grads.d_rv1[3 * j + k] += d_rv1_out[k];
                    mesh_grads.d_rv2[3 * j + k] += d_rv2_out[k];
                }
            }
        }

        // Step 5: Chain-rule MeshGradients -> MeshNodeGradients per surface
        // On first objective, initialize the result structures
        if (obj == 0) {
            result.resize(mesh_input.surfaces.size());
            for (size_t s = 0; s < mesh_input.surfaces.size(); ++s) {
                const auto& sm = mesh_input.surfaces[s];
                result[s].Init(sm.name, sm.nx, sm.ny);
            }
        }

        // Build mapping from input surface index to solver surface index.
        // Each input surface creates one solver surface, plus one mirror if
        // y-duplicated.
        Index solver_surf_idx = 0;
        for (size_t s = 0; s < mesh_input.surfaces.size(); ++s) {
            const auto& sm = mesh_input.surfaces[s];
            const SurfaceMeshInfo& sinfo = mesh.GetSurfaceInfo(solver_surf_idx);

            MeshNodeGradients* target = nullptr;
            switch (obj) {
                case 0:
                    target = &result[s].dCL;
                    break;
                case 1:
                    target = &result[s].dCD;
                    break;
                case 2:
                    target = &result[s].dCY;
                    break;
                case 3:
                    target = &result[s].dCl;
                    break;
                case 4:
                    target = &result[s].dCm;
                    break;
                case 5:
                    target = &result[s].dCn;
                    break;
                default:
                    break;
            }

            // Propagate original surface vortex gradients to nodes
            PropagateMeshGradientsToNodes(mesh_grads, sm,
                                          sinfo.first_vortex_index, sm.nx,
                                          sm.ny, *target);

            ++solver_surf_idx;

            // For y-duplicated surfaces, also propagate the mirror vortex
            // gradients back to the original nodes (with y-reflection and
            // left/right swap)
            if (sm.y_duplicated) {
                const SurfaceMeshInfo& minfo =
                    mesh.GetSurfaceInfo(solver_surf_idx);
                PropagateMirrorMeshGradientsToNodes(
                    mesh_grads, sm, minfo.first_vortex_index,
                    sm.y_duplicate_plane, sm.nx, sm.ny, *target);
                ++solver_surf_idx;
            }
        }
    }

    return result;
}

}  // namespace vlm::ad
