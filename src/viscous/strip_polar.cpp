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

#include "pegasus/viscous/strip_polar.h"

#include <utility>

#include <scirocco/scirocco.h>

namespace pegasus::viscous {

StripPolarSolver::StripPolarSolver(const std::string& naca, int num_panels)
    : solver_(std::make_unique<scirocco::Scirocco>(naca, num_panels)) {}

StripPolarSolver::StripPolarSolver(const Eigen::MatrixXd& coords,
                                   int num_panels)
    : solver_(std::make_unique<scirocco::Scirocco>(coords, num_panels)) {}

StripPolarSolver::~StripPolarSolver() = default;
StripPolarSolver::StripPolarSolver(StripPolarSolver&&) noexcept = default;
StripPolarSolver& StripPolarSolver::operator=(StripPolarSolver&&) noexcept =
    default;

bool StripPolarSolver::SolveFixedAlpha(Scalar alpha_rad, Scalar reynolds,
                                       Scalar mach) {
    const double alpha_deg = static_cast<double>(alpha_rad * kRadiansToDegrees);
    solver_->SetOper(alpha_deg, /*cl=*/-1.0, static_cast<double>(reynolds),
                     static_cast<double>(mach));

    // Warm-start: skip BL reinit if we have a previous converged state for
    // this airfoil.
    solver_->GetState().oper.initbl = !warm_started_;

    const bool converged = solver_->Solve() && solver_->IsConverged();

    if (converged) {
        warm_started_ = true;
        solver_->ComputeSensitivities();
    } else {
        // Failed solve — drop the cached state so the next call retries with
        // a fresh BL initialization.
        warm_started_ = false;
    }
    return converged;
}

bool StripPolarSolver::SolveFixedCl(Scalar cl_target, Scalar reynolds,
                                    Scalar mach) {
    // Scirocco interprets cl ≥ 0 as a target-Cl request and internally
    // solves for the matching α. Seed with the thin-airfoil estimate
    // (α ≈ Cl/(2π)) so Scirocco's Newton starts close to the converged α.
    // For negative target Cl Scirocco's fixed-Cl mode doesn't engage
    // (it requires cl ≥ 0); we punt to fixed-α at the thin-airfoil
    // estimate, which is the right behavior for symmetric airfoils.
    Scalar alpha_seed_deg =
        static_cast<double>(cl_target / (2.0 * kPi)) * kRadiansToDegrees;
    if (cl_target < Scalar(0.0)) {
        return SolveFixedAlpha(cl_target / (2.0 * kPi), reynolds, mach);
    }
    solver_->SetOper(alpha_seed_deg, static_cast<double>(cl_target),
                     static_cast<double>(reynolds), static_cast<double>(mach));
    solver_->GetState().oper.initbl = !warm_started_;
    const bool converged = solver_->Solve() && solver_->IsConverged();
    if (converged) {
        warm_started_ = true;
        solver_->ComputeSensitivities();
    } else {
        warm_started_ = false;
    }
    return converged;
}

Scalar StripPolarSolver::GetAlphaRadians() const {
    return static_cast<Scalar>(solver_->GetAlpha()) * kDegreesToRadians;
}

bool StripPolarSolver::IsConverged() const {
    return solver_->IsConverged();
}

Scalar StripPolarSolver::GetCl() const {
    return static_cast<Scalar>(solver_->GetCl());
}
Scalar StripPolarSolver::GetCd() const {
    return static_cast<Scalar>(solver_->GetCd());
}
Scalar StripPolarSolver::GetCm() const {
    return static_cast<Scalar>(solver_->GetCm());
}

Scalar StripPolarSolver::GetDclDalphaPerRadian() const {
    // Scirocco returns dCl / d(α°); convert to per-radian via chain rule.
    return static_cast<Scalar>(solver_->GetDclDalpha()) * kRadiansToDegrees;
}

Scalar StripPolarSolver::GetDcdDalphaPerRadian() const {
    return static_cast<Scalar>(solver_->GetDcdDalpha()) * kRadiansToDegrees;
}

const Eigen::MatrixXd& StripPolarSolver::GetDclDx() const {
    return solver_->GetDclDx();
}

const Eigen::MatrixXd& StripPolarSolver::GetDcdDx() const {
    return solver_->GetDcdDx();
}

int StripPolarSolver::GetNumPanels() const {
    return static_cast<int>(solver_->GetFoilCoords().cols());
}

}  // namespace pegasus::viscous
