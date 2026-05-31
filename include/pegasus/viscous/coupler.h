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
#include "pegasus/viscous/strip_polar.h"

#include <string>
#include <vector>

namespace pegasus::viscous {

/// @brief Coupler options.
struct CouplerOptions {
    /// NACA 4-digit airfoil applied to every strip (v1 uses a single airfoil
    /// for the whole wing — extend in future versions).
    std::string airfoil_naca = "0012";

    /// Reynolds number used at every strip.
    Scalar reynolds = 1.0e6;

    /// Freestream Mach (passed to Scirocco's compressibility correction).
    Scalar mach = 0.0;

    /// Number of panels in each Scirocco instance.
    int scirocco_panels = 199;

    /// Maximum Newton outer iterations.
    int max_iterations = 20;

    /// Convergence tolerance: ‖R‖∞ in Cl units.
    Scalar tolerance = 1.0e-5;

    /// Newton step damping (line-search factor cap). 1.0 = full step.
    Scalar damping = 1.0;

    /// Maximum |Δα|∞ change per Newton iteration, radians. Step is
    /// uniformly scaled down if exceeded. Prevents overshoot into
    /// non-physical fixed points when J is poorly conditioned (which
    /// happens when the 3D lift slope is close to 2π, making J ≈ small).
    Scalar max_step_rad = 0.5 * kDegreesToRadians;
};

/// @brief Outcome of a coupler solve.
struct CouplerResult {
    bool converged = false;
    int iterations = 0;
    Scalar final_residual = 0.0;
    AerodynamicResults results;
    VecX delta_alpha;  ///< Per-strip Δα at exit, radians.
};

/// @brief Total derivatives of integrated outputs w.r.t. design variables,
/// accounting for the implicit coupling between VLM and 2D viscous solvers.
struct CoupledGradients {
    Scalar CL = 0.0;
    Scalar CD = 0.0;

    /// dCL/dα∞ and dCD/dα∞ (per radian) including the implicit coupling
    /// du/dα∞ that satisfies the converged residual.
    Scalar dCL_dalpha = 0.0;
    Scalar dCD_dalpha = 0.0;
};

/// @brief Viscous-inviscid coupler: drives a VLM ↔ 2D airfoil Newton loop.
///
/// At each outer iteration the coupler:
///   1. Applies the current Δα to every strip's twist_angle
///   2. Re-builds the AIC and re-solves the VLM
///   3. Extracts strip Cl_VLM and α_eff
///   4. Runs Scirocco at α_eff for each strip → Cl_2D
///   5. Residual r_i = Cl_VLM,i − Cl_2D,i
///   6. Newton step with a finite-difference Jacobian column build
///
/// The coupler mutates strip.twist_angle of the supplied solver's mesh.
/// At exit, strip twists are left at the converged state. The original
/// twists can be restored by calling GetBaselineTwists() externally and
/// re-applying them.
class ViscousVlmCoupler {
public:
    ViscousVlmCoupler(VortexLatticeSolver& solver,
                      const FlightCondition& flight,
                      const CouplerOptions& options = {});

    /// @brief Solve the coupled problem.
    CouplerResult Solve();

    /// @brief Total derivatives of CL and CD with respect to the freestream
    /// angle of attack α∞ at the currently-converged state.
    ///
    /// Must be called after `Solve()` has converged. Uses the implicit
    /// function theorem:
    ///   J · (du/dα) = −∂R/∂α     [forward sensitivity]
    ///   dF/dα = ∂F/∂α + ∂F/∂u · (du/dα)
    /// All partials and the Jacobian J = ∂R/∂u are obtained via finite
    /// differences in this v1 implementation.
    ///
    /// On return, the coupler's internal state is left at the converged
    /// (u*, x) baseline.
    CoupledGradients ComputeAlphaDerivatives();

    /// @brief Baseline (pre-coupling) twist of each strip, radians.
    [[nodiscard]] const VecX& GetBaselineTwists() const {
        return baseline_twist_;
    }

    /// @brief Current Δα per strip, radians.
    [[nodiscard]] const VecX& GetDeltaAlpha() const { return delta_alpha_; }

    /// @brief Number of coupling unknowns (= number of strips).
    [[nodiscard]] Index GetNumUnknowns() const {
        return static_cast<Index>(strip_polars_.size());
    }

private:
    /// Apply delta_alpha_ to strip twists and refresh the mesh + AIC.
    void ApplyTwistsAndRefresh();

    /// Run a VLM solve at the current twist state and return strip Cl and
    /// α_eff vectors. Also caches the SolutionState for analytical
    /// linearization use.
    void EvaluateVlm(VecX& cl_vlm, VecX& alpha_eff);

    /// Compute residual R = Cl_VLM − Cl_2D at the current state.
    /// Also fills cd_2d_ and dcl2d_dalpha_/dcd2d_dalpha_ for the
    /// analytical Jacobian and total-derivative chains.
    VecX ComputeResidual(const VecX& alpha_eff, const VecX& cl_vlm);

    /// Build the analytical Jacobian J = ∂R/∂u using a forward
    /// propagation through the VLM chain (one cheap forward solve per
    /// column reusing the cached AIC LU) plus Scirocco's analytical
    /// dCl_2D/dα. No additional VLM or Scirocco solves needed.
    MatX BuildJacobianAnalytical();

    /// Aggregated forward-mode response of a single perturbation
    /// direction through the VLM chain (and Scirocco's linearization).
    struct ForwardResponse {
        VecX dCl_strip;    ///< Per-strip δCl_VLM
        VecX dAlphaEff;    ///< Per-strip δα_eff (induced + explicit)
        Scalar dCL_total;  ///< Integrated δCL (stability axis)
        Scalar dCD_total;  ///< Integrated δCD = δCDi + δCDv
    };

    /// Forward-mode propagation: given a normal-vector perturbation `dn`
    /// (3*N_vortex), a freestream-direction perturbation `dv_inf`, and a
    /// per-strip explicit twist perturbation `dtwist_explicit`, computes
    /// per-strip Cl_VLM / α_eff sensitivities, integrated CL/CD
    /// sensitivities, and the viscous-drag sensitivity through Scirocco.
    ForwardResponse RunForwardChain(const VecX& dn, const Vec3& dv_inf,
                                    const VecX& dtwist_explicit);

    VortexLatticeSolver& solver_;
    FlightCondition flight_;
    CouplerOptions options_;
    ForceCalculator force_calc_;

    /// One 2D solver per VLM strip.
    std::vector<StripPolarSolver> strip_polars_;

    /// Baseline strip.twist_angle captured before any coupling injection.
    VecX baseline_twist_;

    /// Current per-strip Δα (radians).
    VecX delta_alpha_;

    /// Cached 2D drag per strip from the last residual evaluation.
    VecX cd_2d_;
    /// Cached Scirocco dCl_2D/dα per strip (per radian).
    VecX dcl2d_dalpha_;
    /// Cached Scirocco dCd_2D/dα per strip (per radian).
    VecX dcd2d_dalpha_;

    /// Cached SolutionState from the last EvaluateVlm() call; used by the
    /// analytical forward chain. Not valid before the first solve.
    SolutionState last_solution_;

    /// Cached VLM aerodynamic results from the last evaluation.
    AerodynamicResults last_results_;
};

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_COUPLER_H_
