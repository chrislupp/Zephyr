// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/aero/force_calculator.h"

#include <cmath>

namespace vlm {

namespace {

// 1/(2*pi) for 2D Biot-Savart in Trefftz plane
constexpr Scalar kHalfPiInv = 0.15915494309189533577;

}  // namespace

void ForceCalculator::BodyToStabilityAxis(const Vec3& cf_body, Scalar alpha,
                                          Scalar& CL, Scalar& CD) const {
    Scalar ca = std::cos(alpha);
    Scalar sa = std::sin(alpha);

    // Coordinate system: X forward, Y right, Z up
    // Freestream direction at angle alpha: (cos(alpha), 0, sin(alpha))
    // Lift direction (perpendicular, upward): (-sin(alpha), 0, cos(alpha))
    // Drag direction (along freestream): (cos(alpha), 0, sin(alpha))

    // Lift: perpendicular to freestream (positive up)
    // CL = F · L_dir = CZ*cos(alpha) - CX*sin(alpha)
    CL = cf_body.z() * ca - cf_body.x() * sa;

    // Drag: along freestream (positive aft)
    // CD = F · D_dir = CX*cos(alpha) + CZ*sin(alpha)
    CD = cf_body.x() * ca + cf_body.z() * sa;
}

StripLoads ForceCalculator::ComputeStripLoads(
    const VortexLattice& mesh, const SolutionState& solution,
    const FlightCondition& condition, Index strip_index,
    [[maybe_unused]] Scalar reference_area) const {
    StripLoads loads;
    loads.strip_index = strip_index;

    const Strip& strip = mesh.GetStrip(strip_index);

    loads.chord = strip.chord_center;
    loads.width = strip.width;
    loads.area = strip.chord_center * strip.width;
    loads.leading_edge_position = strip.leading_edge_center;
    loads.quarter_chord_position = strip.leading_edge_center;
    loads.quarter_chord_position.x() += 0.25 * strip.chord_center;

    // Freestream velocity direction
    Vec3 v_inf = VortexLatticeSolver::ComputeFreestreamDirection(condition);

    // Accumulate forces from all vortices in strip
    Vec3 total_force = Vec3::Zero();
    Vec3 total_moment = Vec3::Zero();
    Scalar cnc = 0.0;

    for (Index k = 0; k < strip.num_vortices; ++k) {
        Index i = strip.first_vortex_index + k;
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        // Effective velocity at bound vortex = freestream + induced
        Vec3 v_eff = v_inf;
        if (include_induced_velocity_ && solution.induced_velocity.rows() > i) {
            v_eff.x() += solution.induced_velocity(i, 0);
            v_eff.y() += solution.induced_velocity(i, 1);
            v_eff.z() += solution.induced_velocity(i, 2);
        }

        // Bound vortex vector (from left to right)
        Vec3 gamma_vec = vor.bound_vortex_right - vor.bound_vortex_left;

        // Kutta-Joukowski: dF = rho * Gamma * (V x dl)
        // For unit circulation element: F = V_eff x gamma_vec
        // For coefficient: F/(q*S) = 2 * Gamma * (V_eff x gamma_vec) / S
        Scalar gamma = solution.circulation(i);
        Vec3 f = 2.0 * gamma * v_eff.cross(gamma_vec);

        total_force += f;

        // Moment about strip quarter chord
        Vec3 r = vor.control_point - loads.quarter_chord_position;
        Vec3 m = r.cross(f);
        total_moment += m;

        // Span loading contribution (normal force on strip)
        // Using strip normal to project force
        Vec3 strip_normal(0.0, strip.trefftz_normal_y, strip.trefftz_normal_z);
        cnc += strip.chord_center * f.dot(strip_normal);
    }

    // Normalize forces by strip area
    loads.cf = total_force / loads.area;

    // Normalize moments by strip chord * area
    // Convert to aircraft convention (negate x and z for roll/yaw)
    Vec3 cm_aircraft(-total_moment.x(), total_moment.y(), -total_moment.z());
    loads.cm = cm_aircraft / (loads.chord * loads.area);

    // Span loading
    loads.cnc = cnc / loads.area;

    // Compute lift and drag coefficients
    BodyToStabilityAxis(loads.cf, condition.alpha, loads.cl, loads.cd);
    loads.cy = loads.cf.y();

    return loads;
}

AerodynamicResults ForceCalculator::ComputeForces(
    const Aircraft& aircraft, const VortexLattice& mesh,
    const SolutionState& solution, const FlightCondition& condition) const {
    ReferenceQuantities ref;
    ref.reference_area = aircraft.GetReferenceArea();
    ref.reference_chord = aircraft.GetReferenceChord();
    ref.reference_span = aircraft.GetReferenceSpan();
    ref.moment_reference_point = aircraft.GetMomentReferencePoint();
    return ComputeForces(ref, mesh, solution, condition);
}

AerodynamicResults ForceCalculator::ComputeForces(
    const ReferenceQuantities& ref, const VortexLattice& mesh,
    const SolutionState& solution, const FlightCondition& condition) const {
    AerodynamicResults results;

    if (!solution.is_valid || mesh.GetNumVortices() == 0) {
        return results;
    }

    // Set reference quantities
    results.reference_area = ref.reference_area;
    results.reference_chord = ref.reference_chord;
    results.reference_span = ref.reference_span;
    results.moment_reference_point = ref.moment_reference_point;

    // Freestream direction for stability axis transformation (used in
    // BodyToStabilityAxis)

    // Initialize totals
    Vec3 total_force = Vec3::Zero();
    Vec3 total_moment = Vec3::Zero();

    // Process each surface
    results.surface_loads.resize(mesh.GetNumSurfaces());

    for (Index surf_idx = 0; surf_idx < mesh.GetNumSurfaces(); ++surf_idx) {
        const auto& surf_info = mesh.GetSurfaceInfo(surf_idx);
        SurfaceLoads& surf_loads = results.surface_loads[surf_idx];

        surf_loads.surface_index = surf_idx;
        surf_loads.name = surf_info.name;
        surf_loads.area = 0.0;

        Vec3 surf_force = Vec3::Zero();
        Vec3 surf_moment = Vec3::Zero();

        // Find strips belonging to this surface
        for (Index j = 0; j < mesh.GetNumStrips(); ++j) {
            const Strip& strip = mesh.GetStrip(j);
            if (strip.surface_index != surf_idx)
                continue;

            StripLoads strip_loads = ComputeStripLoads(
                mesh, solution, condition, j, results.reference_area);
            surf_loads.strip_loads.push_back(strip_loads);
            surf_loads.area += strip_loads.area;

            // Accumulate surface forces (in strip coordinates, normalized by
            // strip area)
            Vec3 strip_force = strip_loads.cf * strip_loads.area;
            surf_force += strip_force;

            // Moment about surface reference (quarter chord of first strip)
            Vec3 r = strip_loads.quarter_chord_position -
                     results.moment_reference_point;
            surf_moment +=
                strip_loads.cm * strip_loads.chord * strip_loads.area;
            surf_moment += r.cross(strip_force);
        }

        // Normalize surface forces by reference area
        surf_loads.CF = surf_force / results.reference_area;

        // Convert moment to aircraft convention (negate x and z for roll/yaw)
        Vec3 cm_aircraft(-surf_moment.x(), surf_moment.y(), -surf_moment.z());
        surf_loads.CM =
            cm_aircraft / (results.reference_area * results.reference_chord);

        // Convert to lift/drag
        BodyToStabilityAxis(surf_loads.CF, condition.alpha, surf_loads.CL,
                            surf_loads.CD);
        surf_loads.CY = surf_loads.CF.y();

        // Accumulate to totals
        total_force += surf_force;
        total_moment += surf_moment;
    }

    // Total body-axis coefficients
    results.CF = total_force / results.reference_area;
    results.CX = results.CF.x();
    results.CY = results.CF.y();
    results.CZ = results.CF.z();

    // Total moments about reference point
    results.CM =
        total_moment / (results.reference_area * results.reference_chord);

    // Convert from math convention (X fwd, Y right, Z up, right-hand moments)
    // to aircraft convention (positive roll = right wing down, positive yaw =
    // nose right) Pitch (about Y axis) has same sign in both conventions
    results.Cl = -results.CM.x();  // Negate: RH x-moment is left roll, we want
                                   // right roll positive
    results.Cm = results.CM.y();   // Same sign: nose up is positive in both
    results.Cn = -results.CM.z();  // Negate: RH z-moment is nose left, we want
                                   // nose right positive

    // Stability axis coefficients
    BodyToStabilityAxis(results.CF, condition.alpha, results.CL, results.CD);
    results.CDi_nearfield =
        results.CD;  // Near-field drag is all induced for inviscid VLM

    // Compute far-field drag using Trefftz plane method
    ComputeTrefftzPlaneDrag(
        mesh, solution, results.reference_area, results.reference_span,
        results.CL_ff, results.CY_ff, results.CDi_ff, results.span_efficiency);

    results.is_valid = true;
    return results;
}

void ForceCalculator::ComputeTrefftzPlaneDrag(
    const VortexLattice& mesh, const SolutionState& solution,
    Scalar reference_area, Scalar reference_span, Scalar& CL_ff, Scalar& CY_ff,
    Scalar& CDi_ff, Scalar& span_efficiency) const {
    CL_ff = 0.0;
    CY_ff = 0.0;
    CDi_ff = 0.0;
    span_efficiency = 0.0;

    if (!solution.is_valid || mesh.GetNumStrips() == 0) {
        return;
    }

    Index num_strips = mesh.GetNumStrips();

    // Compute total circulation for each strip
    std::vector<Scalar> strip_gamma(num_strips, 0.0);
    std::vector<Vec3> strip_te_left(num_strips);
    std::vector<Vec3> strip_te_right(num_strips);
    std::vector<Vec3> strip_center(num_strips);

    for (Index j = 0; j < num_strips; ++j) {
        const Strip& strip = mesh.GetStrip(j);

        // Sum circulation from all vortices in strip
        for (Index k = 0; k < strip.num_vortices; ++k) {
            Index i = strip.first_vortex_index + k;
            strip_gamma[j] += solution.circulation(i);
        }

        // Get trailing edge vortex (last chordwise vortex)
        Index last_vor = strip.first_vortex_index + strip.num_vortices - 1;
        const HorseshoeVortex& vor = mesh.GetVortex(last_vor);

        // Wake endpoints (trailing vortex legs extend from bound vortex ends)
        strip_te_left[j] = vor.bound_vortex_left;
        strip_te_right[j] = vor.bound_vortex_right;
        strip_center[j] =
            0.5 * (vor.bound_vortex_left + vor.bound_vortex_right);
    }

    // For each strip, compute induced velocity from all other strips
    // in the Trefftz plane (2D Biot-Savart)
    for (Index jc = 0; jc < num_strips; ++jc) {
        const Strip& strip_c = mesh.GetStrip(jc);

        // Strip geometry in Trefftz plane (Y-Z coordinates)
        Scalar dy_c = strip_te_right[jc].y() - strip_te_left[jc].y();
        Scalar dz_c = strip_te_right[jc].z() - strip_te_left[jc].z();
        Scalar ds_c = std::sqrt(dy_c * dy_c + dz_c * dz_c);

        if (ds_c < kEpsilon)
            continue;

        // Normal to strip in Trefftz plane (used for downwash calculation)
        [[maybe_unused]] Scalar ny = -dz_c / ds_c;
        [[maybe_unused]] Scalar nz = dy_c / ds_c;

        // Control point location
        Scalar yc = strip_center[jc].y();
        Scalar zc = strip_center[jc].z();

        // Induced velocity from all strips
        Scalar vy = 0.0;
        Scalar vz = 0.0;

        for (Index jv = 0; jv < num_strips; ++jv) {
            const Strip& strip_v = mesh.GetStrip(jv);

            // Determine core radius
            Scalar ds_v = std::sqrt(
                std::pow(strip_te_right[jv].y() - strip_te_left[jv].y(), 2) +
                std::pow(strip_te_right[jv].z() - strip_te_left[jv].z(), 2));

            // Use surface index to determine if same component
            // (could also use vortex component_index for finer control)
            Scalar core_radius = 0.0;
            if (strip_c.surface_index != strip_v.surface_index) {
                core_radius = std::max(
                    trefftz_core_chord_fraction_ * strip_v.chord_center,
                    trefftz_core_width_fraction_ * ds_v);
            }
            Scalar core_sq = core_radius * core_radius;
            Scalar core_4 = core_sq * core_sq;

            // Distance from control point to strip endpoints
            Scalar dy1 = yc - strip_te_left[jv].y();
            Scalar dy2 = yc - strip_te_right[jv].y();
            Scalar dz1 = zc - strip_te_left[jv].z();
            Scalar dz2 = zc - strip_te_right[jv].z();

            // 2D Biot-Savart with core regularization
            Scalar rsq1 =
                std::sqrt(std::pow(dy1 * dy1 + dz1 * dz1, 2) + core_4);
            Scalar rsq2 =
                std::sqrt(std::pow(dy2 * dy2 + dz2 * dz2, 2) + core_4);

            if (rsq1 > kEpsilon && rsq2 > kEpsilon) {
                vy += kHalfPiInv * strip_gamma[jv] * (dz1 / rsq1 - dz2 / rsq2);
                vz += kHalfPiInv * strip_gamma[jv] * (-dy1 / rsq1 + dy2 / rsq2);
            }
        }

        // Far-field lift (Kutta-Joukowski in Trefftz plane)
        // L = rho * Gamma * span_component_y
        CL_ff += 2.0 * strip_gamma[jc] * dy_c / reference_area;

        // Far-field side force
        CY_ff -= 2.0 * strip_gamma[jc] * dz_c / reference_area;

        // Trefftz plane drag = kinetic energy in crossflow
        // CD = Gamma * (dz*vy - dy*vz) / S_ref
        CDi_ff += strip_gamma[jc] * (dz_c * vy - dy_c * vz) / reference_area;
    }

    // Span efficiency factor
    // e = CL^2 / (pi * AR * CDi)
    Scalar AR = reference_span * reference_span / reference_area;
    if (std::abs(CDi_ff) > kEpsilon) {
        span_efficiency = (CL_ff * CL_ff + CY_ff * CY_ff) / (kPi * AR * CDi_ff);
    }
}

}  // namespace vlm
