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

#ifndef PEGASUS_VISCOUS_COUPLER_H_
#define PEGASUS_VISCOUS_COUPLER_H_

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/aero/force_calculator.h"
#include "pegasus/core/types.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/viscous/coupling_station.h"
#include "pegasus/viscous/strip_polar.h"

#include <vector>

namespace pegasus::viscous {

/// @brief Coupler options.
///
/// `stations` defines the user-controlled 2D viscous evaluation grid. The
/// VLM mesh is left untouched (the coupler only mutates strip twist
/// offsets). The number, spanwise placement (eta), airfoil shape, panel
/// count, and Reynolds number are all per-station.
struct CouplerOptions {
    /// User-supplied 2D viscous stations (≥ 2, strictly increasing eta).
    std::vector<CouplingStation> stations;

    /// Freestream Mach (passed to Scirocco's compressibility correction).
    Scalar mach = 0.0;

    /// Maximum Newton outer iterations.
    int max_iterations = 20;

    /// Convergence tolerance: ‖R‖∞ in Cl units.
    Scalar tolerance = 1.0e-5;

    /// Newton step damping (line-search factor cap). 1.0 = full step.
    Scalar damping = 1.0;

    /// Maximum |Δα|∞ change per Newton iteration, radians.
    Scalar max_step_rad = 1.0 * kDegreesToRadians;
};

/// @brief Outcome of a coupler solve.
struct CouplerResult {
    bool converged = false;
    int iterations = 0;
    Scalar final_residual = 0.0;
    AerodynamicResults results;

    /// Per-station Δα at exit, radians.
    VecX delta_alpha_station;

    /// Per-strip twist offset (= spline_to_strip · delta_alpha_station),
    /// radians. This is what was actually written to mesh strip twists.
    VecX twist_offset_strip;
};

/// @brief Total derivatives of integrated outputs w.r.t. design variables,
/// accounting for the implicit coupling.
struct CoupledGradients {
    Scalar CL = 0.0;
    Scalar CD = 0.0;

    /// dCL/dα∞ and dCD/dα∞ (per radian) including the implicit coupling
    /// du/dα∞ that satisfies the converged residual.
    Scalar dCL_dalpha = 0.0;
    Scalar dCD_dalpha = 0.0;
};

/// @brief Viscous-inviscid coupler: drives a VLM ↔ 2D airfoil Newton loop
/// using a user-defined station grid that is independent of the VLM
/// strips.
///
/// Workflow per outer iteration:
///   1. Δα_station ⟶ twist_offset_strip via cubic-spline `to_strip` op
///   2. Strip twists set, mesh normals refreshed, AIC rebuilt, VLM solved
///   3. Per-strip α_eff and Cl_local computed in each strip's intrinsic
///      local frame (chord direction, spanwise tangent, lift normal).
///   4. Spline-interpolate strip α_eff and Cl_local to station values via
///      the `to_station` op.
///   5. Run Scirocco at each station's α_eff (fixed-α mode, station's own
///      airfoil coords, station's Reynolds).
///   6. Residual R_k = Cl_local_station_k − Cl_2D_station_k.
///   7. Newton step (current v: FD Jacobian; analytical to come).
///
/// Sweep / dihedral / per-station airfoils all fall out of this structure
/// because every projection happens in the strip's actual 3D local frame,
/// and the station grid is purely user-defined eta values.
class ViscousVlmCoupler {
public:
    ViscousVlmCoupler(VortexLatticeSolver& solver,
                      const FlightCondition& flight,
                      const CouplerOptions& options);

    /// @brief Solve the coupled problem.
    CouplerResult Solve();

    /// @brief Total derivatives of CL and CD with respect to the freestream
    /// angle of attack α∞ at the currently-converged state. Must be called
    /// after `Solve()` has converged.
    CoupledGradients ComputeAlphaDerivatives();

    /// @brief Baseline (pre-coupling) twist of each strip, radians.
    [[nodiscard]] const VecX& GetBaselineTwists() const {
        return baseline_twist_;
    }

    /// @brief Current per-station Δα, radians.
    [[nodiscard]] const VecX& GetDeltaAlpha() const { return delta_alpha_; }

    /// @brief Number of coupling unknowns (= number of stations).
    [[nodiscard]] Index GetNumUnknowns() const {
        return static_cast<Index>(station_polars_.size());
    }

    /// @brief Strip↔station spline operators built at construction.
    [[nodiscard]] const StationInterpolation& GetInterpolation() const {
        return interp_;
    }

private:
    /// Apply current per-station Δα to strip twists via the spline operator
    /// and refresh the mesh + AIC.
    void ApplyTwistsAndRefresh();

    /// Run a VLM solve at the current twist state. Caches the resulting
    /// SolutionState and AerodynamicResults; returns per-strip α_eff and
    /// Cl_local computed in each strip's intrinsic local frame.
    void EvaluateVlm(VecX& cl_local_strip, VecX& alpha_eff_strip);

    /// Compute residual R = Cl_local_station − Cl_2D(α_eff_station) by
    /// interpolating per-strip quantities to stations and calling each
    /// station's Scirocco solver. Caches Scirocco linearization values.
    VecX ComputeResidual(const VecX& alpha_eff_strip,
                         const VecX& cl_local_strip);

    /// Aggregated forward-mode response through the VLM chain + Scirocco
    /// linearization for a single perturbation direction.
    struct ForwardResponse {
        VecX dCl_station;        ///< δCl_VLM_local per station
        VecX dAlphaEff_station;  ///< δα_eff per station (incl. −2·twist term)
        Scalar dCL_total;        ///< δCL_total (stability axis)
        Scalar dCD_total;        ///< δCD_total (CDi_nearfield + CDv)
    };

    /// Forward chain. Inputs:
    ///   dtwist_strip — per-strip twist perturbation (radians). For a
    ///                  station-Δα perturbation this is `to_strip · e_k`;
    ///                  zero for an α∞ perturbation.
    ///   dv_inf       — freestream-direction perturbation. Zero for a
    ///                  station-Δα perturbation; df_inf/dα · dα for α∞.
    ForwardResponse RunForwardChain(const VecX& dtwist_strip,
                                    const Vec3& dv_inf);

    /// Build the analytical Newton Jacobian J = ∂R/∂(Δα_station) by
    /// running one forward chain per station.
    MatX BuildJacobianAnalytical();

    VortexLatticeSolver& solver_;
    FlightCondition flight_;
    CouplerOptions options_;
    ForceCalculator force_calc_;
    StationInterpolation interp_;

    /// One 2D solver per coupling station (not per strip).
    std::vector<StripPolarSolver> station_polars_;

    /// Baseline strip.twist_angle captured before any coupling injection,
    /// radians.
    VecX baseline_twist_;

    /// Current per-station Δα, radians.
    VecX delta_alpha_;

    /// Cached per-station 2D drag from the last residual evaluation.
    VecX cd_2d_station_;
    /// Cached per-station Scirocco dCl_2D/dα (per radian).
    VecX dcl2d_dalpha_station_;
    /// Cached per-station Scirocco dCd_2D/dα (per radian).
    VecX dcd2d_dalpha_station_;

    /// Cached SolutionState from the last EvaluateVlm() call.
    SolutionState last_solution_;
    /// Cached VLM aerodynamic results from the last evaluation.
    AerodynamicResults last_results_;
};

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_COUPLER_H_
