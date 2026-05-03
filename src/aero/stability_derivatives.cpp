/*
    Project Pegasus
    Version 0.1.0

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

#include "pegasus/aero/stability_derivatives.h"

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/aero/force_calculator.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <cmath>

namespace pegasus {

namespace {

/// Compute forces from a unit circulation distribution
AerodynamicResults ComputeUnitForces(const Aircraft& aircraft,
                                     const VortexLattice& mesh,
                                     const VecX& unit_circulation,
                                     const FlightCondition& condition) {
    ForceCalculator calculator;

    // Build a solution state with the unit circulation
    SolutionState unit_solution;
    unit_solution.circulation = unit_circulation;
    unit_solution.is_valid = true;

    // Don't include induced velocity for unit response (already in base
    // solution)
    calculator.SetIncludeInducedVelocity(false);

    return calculator.ComputeForces(aircraft, mesh, unit_solution, condition);
}

}  // namespace

StabilityDerivatives StabilityDerivativeCalculator::ComputeDerivatives(
    const Aircraft& aircraft, VortexLatticeSolver& solver,
    const FlightCondition& condition) {
    StabilityDerivatives derivs;

    // Ensure unit responses are computed
    if (!solver.AreUnitResponsesValid()) {
        solver.SolveUnitResponses();
    }

    // Get the solution state which contains unit circulations
    const SolutionState& solution = solver.GetSolution();
    const VortexLattice& mesh = solver.GetMesh();

    // Reference quantities
    derivs.reference_area = aircraft.GetReferenceArea();
    derivs.reference_chord = aircraft.GetReferenceChord();
    derivs.reference_span = aircraft.GetReferenceSpan();
    derivs.velocity = 1.0;  // Normalized
    derivs.alpha = condition.alpha;
    derivs.beta = condition.beta;

    // Compute forces for each unit response
    // These give dC/d(velocity_component)
    AerodynamicResults F_u = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_u, condition);
    AerodynamicResults F_v = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_v, condition);
    AerodynamicResults F_w = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_w, condition);
    AerodynamicResults F_p = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_p, condition);
    AerodynamicResults F_q = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_q, condition);
    AerodynamicResults F_r = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_r, condition);

    // Convert velocity derivatives to angle derivatives
    // V_inf = V * [cos(alpha)cos(beta), sin(beta), sin(alpha)cos(beta)]
    // dV_inf/dalpha = V * [-sin(alpha)cos(beta), 0, cos(alpha)cos(beta)]
    // dV_inf/dbeta  = V * [-cos(alpha)sin(beta), cos(beta),
    // -sin(alpha)sin(beta)]
    //
    // At beta = 0:
    // dV/dalpha = [-sin(alpha), 0, cos(alpha)]
    // dV/dbeta  = [0, 1, 0]

    Scalar ca = std::cos(condition.alpha);
    Scalar sa = std::sin(condition.alpha);
    Scalar cb = std::cos(condition.beta);
    Scalar sb = std::sin(condition.beta);

    // Chain rule: dC/dalpha = dC/du * du/dalpha + dC/dw * dw/dalpha
    // At small beta: dC/dalpha ≈ -sa * dC/du + ca * dC/dw
    Scalar du_da = -sa * cb;
    Scalar dw_da = ca * cb;
    Scalar du_db = -ca * sb;
    Scalar dv_db = cb;
    Scalar dw_db = -sa * sb;

    // Angle of attack derivatives (stability axis)
    derivs.CL_alpha = F_u.CL * du_da + F_w.CL * dw_da;
    derivs.CY_alpha = F_u.CY * du_da + F_w.CY * dw_da;
    derivs.CD_alpha = F_u.CD * du_da + F_w.CD * dw_da;
    derivs.Cl_alpha = F_u.Cl * du_da + F_w.Cl * dw_da;
    derivs.Cm_alpha = F_u.Cm * du_da + F_w.Cm * dw_da;
    derivs.Cn_alpha = F_u.Cn * du_da + F_w.Cn * dw_da;

    // Sideslip derivatives (stability axis)
    derivs.CL_beta = F_u.CL * du_db + F_v.CL * dv_db + F_w.CL * dw_db;
    derivs.CY_beta = F_u.CY * du_db + F_v.CY * dv_db + F_w.CY * dw_db;
    derivs.CD_beta = F_u.CD * du_db + F_v.CD * dv_db + F_w.CD * dw_db;
    derivs.Cl_beta = F_u.Cl * du_db + F_v.Cl * dv_db + F_w.Cl * dw_db;
    derivs.Cm_beta = F_u.Cm * du_db + F_v.Cm * dv_db + F_w.Cm * dw_db;
    derivs.Cn_beta = F_u.Cn * du_db + F_v.Cn * dv_db + F_w.Cn * dw_db;

    // Rate derivatives (already non-dimensionalized in solver)
    // Roll rate derivatives (p * b / 2V)
    derivs.CL_p = F_p.CL;
    derivs.CY_p = F_p.CY;
    derivs.CD_p = F_p.CD;
    derivs.Cl_p = F_p.Cl;
    derivs.Cm_p = F_p.Cm;
    derivs.Cn_p = F_p.Cn;

    // Pitch rate derivatives (q * c / 2V)
    derivs.CL_q = F_q.CL;
    derivs.CY_q = F_q.CY;
    derivs.CD_q = F_q.CD;
    derivs.Cl_q = F_q.Cl;
    derivs.Cm_q = F_q.Cm;
    derivs.Cn_q = F_q.Cn;

    // Yaw rate derivatives (r * b / 2V)
    derivs.CL_r = F_r.CL;
    derivs.CY_r = F_r.CY;
    derivs.CD_r = F_r.CD;
    derivs.Cl_r = F_r.Cl;
    derivs.Cm_r = F_r.Cm;
    derivs.Cn_r = F_r.Cn;

    // Compute neutral point (location where Cm_alpha = 0)
    // x_np = x_ref - Cm_alpha / CL_alpha * c_ref
    if (std::abs(derivs.CL_alpha) > kEpsilon) {
        derivs.neutral_point_x =
            aircraft.GetMomentReferencePoint().x() -
            derivs.Cm_alpha / derivs.CL_alpha * derivs.reference_chord;
        derivs.static_margin =
            (derivs.neutral_point_x - aircraft.GetMomentReferencePoint().x()) /
            derivs.reference_chord;
    }

    derivs.is_valid = true;
    return derivs;
}

BodyAxisDerivatives StabilityDerivativeCalculator::ComputeBodyAxisDerivatives(
    const Aircraft& aircraft, VortexLatticeSolver& solver,
    const FlightCondition& condition) {
    BodyAxisDerivatives derivs;

    // Ensure unit responses are computed
    if (!solver.AreUnitResponsesValid()) {
        solver.SolveUnitResponses();
    }

    const SolutionState& solution = solver.GetSolution();
    const VortexLattice& mesh = solver.GetMesh();

    // Reference quantities
    derivs.reference_area = aircraft.GetReferenceArea();
    derivs.reference_chord = aircraft.GetReferenceChord();
    derivs.reference_span = aircraft.GetReferenceSpan();
    derivs.velocity = 1.0;
    derivs.alpha = condition.alpha;
    derivs.beta = condition.beta;

    // Compute forces for each unit response (body axis coefficients)
    AerodynamicResults F_u = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_u, condition);
    AerodynamicResults F_v = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_v, condition);
    AerodynamicResults F_w = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_w, condition);
    AerodynamicResults F_p = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_p, condition);
    AerodynamicResults F_q = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_q, condition);
    AerodynamicResults F_r = ComputeUnitForces(
        aircraft, mesh, solution.circulation_per_r, condition);

    // Velocity to angle transformation (same as stability axis)
    Scalar ca = std::cos(condition.alpha);
    Scalar sa = std::sin(condition.alpha);
    Scalar cb = std::cos(condition.beta);
    Scalar sb = std::sin(condition.beta);

    Scalar du_da = -sa * cb;
    Scalar dw_da = ca * cb;
    Scalar du_db = -ca * sb;
    Scalar dv_db = cb;
    Scalar dw_db = -sa * sb;

    // Angle of attack derivatives (body axis)
    derivs.CX_alpha = F_u.CX * du_da + F_w.CX * dw_da;
    derivs.CY_alpha = F_u.CY * du_da + F_w.CY * dw_da;
    derivs.CZ_alpha = F_u.CZ * du_da + F_w.CZ * dw_da;
    derivs.Cl_alpha = F_u.Cl * du_da + F_w.Cl * dw_da;
    derivs.Cm_alpha = F_u.Cm * du_da + F_w.Cm * dw_da;
    derivs.Cn_alpha = F_u.Cn * du_da + F_w.Cn * dw_da;

    // Sideslip derivatives (body axis)
    derivs.CX_beta = F_u.CX * du_db + F_v.CX * dv_db + F_w.CX * dw_db;
    derivs.CY_beta = F_u.CY * du_db + F_v.CY * dv_db + F_w.CY * dw_db;
    derivs.CZ_beta = F_u.CZ * du_db + F_v.CZ * dv_db + F_w.CZ * dw_db;
    derivs.Cl_beta = F_u.Cl * du_db + F_v.Cl * dv_db + F_w.Cl * dw_db;
    derivs.Cm_beta = F_u.Cm * du_db + F_v.Cm * dv_db + F_w.Cm * dw_db;
    derivs.Cn_beta = F_u.Cn * du_db + F_v.Cn * dv_db + F_w.Cn * dw_db;

    // Rate derivatives (body axis)
    derivs.CX_p = F_p.CX;
    derivs.CY_p = F_p.CY;
    derivs.CZ_p = F_p.CZ;
    derivs.Cl_p = F_p.Cl;
    derivs.Cm_p = F_p.Cm;
    derivs.Cn_p = F_p.Cn;

    derivs.CX_q = F_q.CX;
    derivs.CY_q = F_q.CY;
    derivs.CZ_q = F_q.CZ;
    derivs.Cl_q = F_q.Cl;
    derivs.Cm_q = F_q.Cm;
    derivs.Cn_q = F_q.Cn;

    derivs.CX_r = F_r.CX;
    derivs.CY_r = F_r.CY;
    derivs.CZ_r = F_r.CZ;
    derivs.Cl_r = F_r.Cl;
    derivs.Cm_r = F_r.Cm;
    derivs.Cn_r = F_r.Cn;

    derivs.is_valid = true;
    return derivs;
}

}  // namespace pegasus
