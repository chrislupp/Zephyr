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

#ifndef PEGASUS_AERO_FORCE_CALCULATOR_H_
#define PEGASUS_AERO_FORCE_CALCULATOR_H_

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/core/types.h"
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

namespace pegasus {

/// Per-vortex force and moment data in body frame.
///
/// Stores the raw (un-normalized) per-vortex forces before any strip or
/// surface aggregation.  The force on vortex i occupies
/// forces.segment<3>(3*i).  These are the raw Kutta-Joukowski values:
///   F_i = 2 * Gamma_i * (V_eff_i x dl_i)
/// which equal the force coefficient times reference area in the Pegasus
/// non-dimensionalization.  Perseus handles dimensional scaling externally.
struct PerVortexForces {
    VecX forces;   ///< (3*N x 1) body-frame force per vortex
    VecX moments;  ///< (3*N x 1) body-frame moment per vortex (about control
                   ///< point; currently zero — moment arm handled by Perseus)
};

/// Computes aerodynamic forces and moments using near-field integration.
///
/// Uses the Kutta-Joukowski theorem to compute forces on each vortex element:
///   F = rho * V_eff × Gamma
///
/// where V_eff is the total velocity at the bound vortex (freestream +
/// rotation + induced).
class ForceCalculator {
public:
    ForceCalculator() = default;

    /// Compute forces from a solved vortex lattice solution.
    ///
    /// @param aircraft The aircraft geometry (for reference quantities)
    /// @param mesh The vortex lattice mesh
    /// @param solution The converged circulation solution
    /// @param condition The flight condition used in the solve
    /// @return Complete aerodynamic results with forces and moments
    [[nodiscard]] AerodynamicResults ComputeForces(
        const Aircraft& aircraft, const VortexLattice& mesh,
        const SolutionState& solution, const FlightCondition& condition) const;

    /// Compute forces using reference quantities directly (mesh-based input).
    [[nodiscard]] AerodynamicResults ComputeForces(
        const ReferenceQuantities& ref, const VortexLattice& mesh,
        const SolutionState& solution, const FlightCondition& condition) const;

    /// Compute far-field (Trefftz plane) induced drag.
    ///
    /// Uses kinetic energy integral in the wake to compute induced drag.
    /// More accurate than near-field integration for induced drag.
    ///
    /// @param mesh The vortex lattice mesh
    /// @param solution The converged circulation solution
    /// @param reference_area Reference wing area for coefficient normalization
    /// @param reference_span Reference span for efficiency calculation
    /// @param[out] CL_ff Far-field lift coefficient
    /// @param[out] CY_ff Far-field side force coefficient
    /// @param[out] CDi_ff Far-field induced drag coefficient
    /// @param[out] span_efficiency Oswald span efficiency factor
    void ComputeTrefftzPlaneDrag(const VortexLattice& mesh,
                                 const SolutionState& solution,
                                 Scalar reference_area, Scalar reference_span,
                                 Scalar& CL_ff, Scalar& CY_ff, Scalar& CDi_ff,
                                 Scalar& span_efficiency) const;

    /// Compute per-vortex forces in body frame without strip aggregation.
    ///
    /// Returns the raw Kutta-Joukowski force for every horseshoe vortex
    /// element before any strip or surface aggregation.  The force on
    /// vortex i is stored in result.forces.segment<3>(3*i).
    ///
    /// Moments are set to zero — the moment arm about a structural node is
    /// handled by the Perseus load-transfer module.
    ///
    /// @param mesh      The vortex lattice mesh
    /// @param solution  The converged circulation solution
    /// @param condition The flight condition used in the solve
    /// @return PerVortexForces with forces (3N x 1) and moments (3N x 1)
    [[nodiscard]] PerVortexForces ComputePerVortexForces(
        const VortexLattice& mesh, const SolutionState& solution,
        const FlightCondition& condition) const;

    /// Set whether to include induced velocity in force computation.
    /// Default is true (include induced velocity for more accurate forces).
    void SetIncludeInducedVelocity(bool include) {
        include_induced_velocity_ = include;
    }

    /// Set core radius fractions for Trefftz plane induced velocity
    /// computation.
    void SetTrefftzCoreRadius(Scalar chord_fraction, Scalar width_fraction) {
        trefftz_core_chord_fraction_ = chord_fraction;
        trefftz_core_width_fraction_ = width_fraction;
    }

private:
    /// Compute strip-level forces and moments.
    [[nodiscard]] StripLoads ComputeStripLoads(const VortexLattice& mesh,
                                               const SolutionState& solution,
                                               const FlightCondition& condition,
                                               Index strip_index,
                                               Scalar reference_area) const;

    /// Transform body-axis forces to stability-axis forces.
    void BodyToStabilityAxis(const Vec3& cf_body, Scalar alpha, Scalar& CL,
                             Scalar& CD) const;

    bool include_induced_velocity_{true};
    Scalar trefftz_core_chord_fraction_{0.01};
    Scalar trefftz_core_width_fraction_{0.1};
};

}  // namespace pegasus

#endif  // PEGASUS_AERO_FORCE_CALCULATOR_H_
