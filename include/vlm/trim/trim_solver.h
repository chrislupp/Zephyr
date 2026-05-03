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

#ifndef VLM_TRIM_TRIM_SOLVER_H_
#define VLM_TRIM_TRIM_SOLVER_H_

#include "vlm/aero/aerodynamic_results.h"
#include "vlm/aero/stability_derivatives.h"
#include "vlm/core/types.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/solver/vortex_lattice_solver.h"
#include "vlm/trim/flight_constraints.h"

#include <string>
#include <vector>

namespace vlm {

/// Results from trim solution
struct TrimResult {
    bool converged{false};      ///< Whether trim converged
    int iterations{0};          ///< Number of iterations used
    Scalar residual_norm{0.0};  ///< Final residual norm

    /// Final operating state
    OperatingVariables state;

    /// Aerodynamic results at trim condition
    AerodynamicResults aero;

    /// Stability derivatives at trim condition
    StabilityDerivatives derivatives;

    /// Error message if failed
    std::string error_message;
};

/// Options for trim solver
struct TrimOptions {
    int max_iterations{20};    ///< Maximum Newton iterations
    Scalar tolerance{1.0e-6};  ///< Convergence tolerance
    Scalar max_alpha{1.5};     ///< Maximum alpha magnitude (rad, ~86 deg)
    Scalar max_beta{1.5};      ///< Maximum beta magnitude (rad, ~86 deg)
    Scalar max_rate{5.0};      ///< Maximum rate magnitude (non-dim)
    Scalar relaxation{1.0};    ///< Under-relaxation factor (0-1)
    bool verbose{false};       ///< Print iteration info
};

/// Newton-Raphson trim solver for VLM
///
/// Adjusts operating variables (alpha, beta, rates, controls) to satisfy
/// specified constraints (CL, moments, etc.) using Newton-Raphson iteration.
///
/// The solver builds a Jacobian matrix relating changes in operating variables
/// to changes in aerodynamic outputs, then iteratively solves:
///   J * delta_x = -residual
///
/// until the residual (difference between computed and target values) is
/// below the specified tolerance.
///
/// Example usage:
/// @code
/// Aircraft aircraft = LoadAircraft("aircraft.avl");
/// VortexLatticeSolver solver;
/// solver.SetAircraft(aircraft);
///
/// // Set up trim for level flight at CL = 0.5
/// TrimConstraints constraints = TrimPresets::LevelFlight(0.5);
///
/// TrimSolver trimmer;
/// TrimResult result = trimmer.Solve(aircraft, solver, constraints);
///
/// if (result.converged) {
///   std::cout << "Trimmed alpha = " << result.state.alpha << " rad\n";
/// }
/// @endcode
class TrimSolver {
public:
    TrimSolver() = default;

    /// Set solver options
    void SetOptions(const TrimOptions& options) { options_ = options; }

    /// Get current options
    [[nodiscard]] const TrimOptions& GetOptions() const { return options_; }

    /// Solve for trim condition
    ///
    /// @param aircraft Aircraft geometry
    /// @param solver VLM solver (will be modified during iteration)
    /// @param constraints Trim constraints specifying targets
    /// @param initial_state Optional initial guess for operating state
    /// @return TrimResult with converged state or error
    TrimResult Solve(
        const Aircraft& aircraft, VortexLatticeSolver& solver,
        const TrimConstraints& constraints,
        const OperatingVariables& initial_state = OperatingVariables());

private:
    /// Build the residual vector for current state
    [[nodiscard]] VecX BuildResidual(const AerodynamicResults& aero,
                                     const OperatingVariables& state,
                                     const TrimConstraints& constraints) const;

    /// Build the Jacobian matrix
    [[nodiscard]] MatX BuildJacobian(const StabilityDerivatives& derivs,
                                     const OperatingVariables& state,
                                     const TrimConstraints& constraints) const;

    /// Apply constraints to build the system
    void BuildSystem(const AerodynamicResults& aero,
                     const StabilityDerivatives& derivs,
                     const OperatingVariables& state,
                     const TrimConstraints& constraints, MatX& jacobian,
                     VecX& residual) const;

    /// Update state with Newton step
    bool UpdateState(OperatingVariables& state, const VecX& delta,
                     const TrimConstraints& constraints) const;

    /// Check if state is within bounds
    [[nodiscard]] bool CheckBounds(const OperatingVariables& state) const;

    /// Get index of variable in state vector
    [[nodiscard]] Index GetVariableIndex(
        const TrimConstraints& constraints, int var_type,
        const std::string& control_name = "") const;

    TrimOptions options_;
};

}  // namespace vlm

#endif  // VLM_TRIM_TRIM_SOLVER_H_
