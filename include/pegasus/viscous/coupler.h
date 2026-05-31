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

    /// Finite-difference step for Jacobian column assembly, in radians.
    Scalar fd_step_rad = 0.05 * kDegreesToRadians;

    /// Newton step damping (line-search factor cap). 1.0 = full step.
    Scalar damping = 1.0;
};

/// @brief Outcome of a coupler solve.
struct CouplerResult {
    bool converged = false;
    int iterations = 0;
    Scalar final_residual = 0.0;
    AerodynamicResults results;
    VecX delta_alpha;  ///< Per-strip Δα at exit, radians.
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
    /// α_eff vectors.
    void EvaluateVlm(VecX& cl_vlm, VecX& alpha_eff);

    /// Compute residual R = Cl_VLM − Cl_2D at the current state.
    /// Also fills cd_2d_ for drag bookkeeping.
    VecX ComputeResidual(const VecX& alpha_eff, const VecX& cl_vlm);

    /// Build the Jacobian by finite differences over delta_alpha_.
    MatX BuildJacobianFD(const VecX& base_residual);

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

    /// Cached VLM aerodynamic results from the last evaluation.
    AerodynamicResults last_results_;
};

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_COUPLER_H_
