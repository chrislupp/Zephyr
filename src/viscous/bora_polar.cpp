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

#include "pegasus/viscous/bora_polar.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <utility>

#include <bora/kulfan.h>
#include <bora/neuralfoil.h>
#include <bora/types.h>

namespace pegasus::viscous {

namespace {

bora::ModelSize ToBoraModelSize(BoraModelSize size) {
    switch (size) {
        case BoraModelSize::kXXSmall:
            return bora::ModelSize::kXXSmall;
        case BoraModelSize::kXSmall:
            return bora::ModelSize::kXSmall;
        case BoraModelSize::kSmall:
            return bora::ModelSize::kSmall;
        case BoraModelSize::kMedium:
            return bora::ModelSize::kMedium;
        case BoraModelSize::kLarge:
            return bora::ModelSize::kLarge;
        case BoraModelSize::kXLarge:
            return bora::ModelSize::kXLarge;
        case BoraModelSize::kXXLarge:
            return bora::ModelSize::kXXLarge;
        case BoraModelSize::kXXXLarge:
            return bora::ModelSize::kXXXLarge;
    }
    return bora::ModelSize::kXLarge;
}

/// Convert a 2×N airfoil node matrix into normalized Selig-order points for
/// Bora's Kulfan fit. Selig order (upper TE → LE → lower TE) traverses the
/// loop counter-clockwise for the standard x-right / y-up convention, so a
/// clockwise input (e.g. Scirocco's) is detected via the shoelace signed
/// area and reversed. The points are then chord-normalized: LE shifted to
/// the origin, chord scaled to 1.
std::vector<bora::Point2> ToSeligNormalized(const Eigen::MatrixXd& coords) {
    if (coords.rows() != 2 || coords.cols() < 5) {
        throw std::invalid_argument(
            "BoraPolarSolver: coords must be a 2xN matrix with N >= 5.");
    }

    Eigen::Index n = coords.cols();
    // Drop a duplicated closing node so the TE appears exactly twice at the
    // ends after any reversal.
    if ((coords.col(0) - coords.col(n - 1)).norm() < 1e-12) {
        --n;
    }

    std::vector<bora::Point2> pts;
    pts.reserve(static_cast<std::size_t>(n) + 1);
    for (Eigen::Index j = 0; j < n; ++j) {
        pts.push_back({coords(0, j), coords(1, j)});
    }

    double signed_area = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const bora::Point2& a = pts[i];
        const bora::Point2& b = pts[(i + 1) % pts.size()];
        signed_area += a[0] * b[1] - b[0] * a[1];
    }
    if (signed_area < 0.0) {
        std::reverse(pts.begin(), pts.end());
    }

    // Chord-normalize: LE (minimum-x node) to the origin, chord to 1.
    std::size_t le = 0;
    double x_min = pts[0][0];
    double x_max = pts[0][0];
    for (std::size_t i = 1; i < pts.size(); ++i) {
        if (pts[i][0] < x_min) {
            x_min = pts[i][0];
            le = i;
        }
        x_max = std::max(x_max, pts[i][0]);
    }
    const double chord = x_max - x_min;
    if (chord < 1e-12) {
        throw std::invalid_argument(
            "BoraPolarSolver: degenerate airfoil coordinates (zero chord).");
    }
    const double y_le = pts[le][1];
    for (bora::Point2& p : pts) {
        p[0] = (p[0] - x_min) / chord;
        p[1] = (p[1] - y_le) / chord;
    }
    return pts;
}

}  // namespace

struct BoraPolarSolver::Impl {
    std::shared_ptr<const bora::NeuralFoil> model;
    BoraOptions options{};
    bora::KulfanParameters kulfan{};
    bora::AeroResult result{};
    Scalar alpha_rad = 0.0;
    bool converged = false;
    Eigen::MatrixXd dcl_dshape = Eigen::MatrixXd::Zero(kBoraNumShapeParams, 1);
    Eigen::MatrixXd dcd_dshape = Eigen::MatrixXd::Zero(kBoraNumShapeParams, 1);

    void Initialize(std::shared_ptr<const bora::NeuralFoil> shared_model,
                    const Eigen::MatrixXd& coords, const BoraOptions& opts) {
        model = std::move(shared_model);
        options = opts;
        const std::vector<bora::Point2> pts = ToSeligNormalized(coords);
        kulfan = bora::FitKulfanParameters(pts);
    }
};

BoraPolarSolver::BoraPolarSolver() : impl_(std::make_unique<Impl>()) {}

BoraPolarSolver::BoraPolarSolver(const Eigen::MatrixXd& coords,
                                 const BoraOptions& options)
    : impl_(std::make_unique<Impl>()) {
    auto model = std::make_shared<const bora::NeuralFoil>(
        ToBoraModelSize(options.model_size),
        std::filesystem::path(options.data_dir));
    impl_->Initialize(std::move(model), coords, options);
}

BoraPolarSolver::~BoraPolarSolver() = default;
BoraPolarSolver::BoraPolarSolver(BoraPolarSolver&&) noexcept = default;
BoraPolarSolver& BoraPolarSolver::operator=(BoraPolarSolver&&) noexcept =
    default;

bool BoraPolarSolver::SolveFixedAlpha(Scalar alpha_rad, Scalar reynolds,
                                      Scalar mach) {
    // The surrogate is trained on incompressible XFoil data and has no Mach
    // input; callers should treat results as low-subsonic section data.
    (void)mach;

    bora::AeroInput input;
    input.kulfan = impl_->kulfan;
    input.alpha = static_cast<double>(alpha_rad * kRadiansToDegrees);
    input.Re = static_cast<double>(reynolds);
    input.n_crit = static_cast<double>(impl_->options.n_crit);
    input.xtr_upper = static_cast<double>(impl_->options.xtr_upper);
    input.xtr_lower = static_cast<double>(impl_->options.xtr_lower);

    impl_->result = impl_->model->AnalyzeWithGradients(input);
    impl_->alpha_rad = alpha_rad;
    impl_->converged = impl_->result.output.analysis_confidence >=
                       impl_->options.min_confidence;

    for (int i = 0; i < kBoraNumShapeParams; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        impl_->dcl_dshape(i, 0) = impl_->result.gradient.CL[idx];
        impl_->dcd_dshape(i, 0) = impl_->result.gradient.CD[idx];
    }
    return impl_->converged;
}

Scalar BoraPolarSolver::GetAlphaRadians() const {
    return impl_->alpha_rad;
}

bool BoraPolarSolver::IsConverged() const {
    return impl_->converged;
}

Scalar BoraPolarSolver::GetCl() const {
    return impl_->result.output.CL;
}

Scalar BoraPolarSolver::GetCd() const {
    return impl_->result.output.CD;
}

Scalar BoraPolarSolver::GetCm() const {
    return impl_->result.output.CM;
}

Scalar BoraPolarSolver::GetDclDalphaPerRadian() const {
    // Bora's α gradient is per-degree (its α input is in degrees).
    return impl_->result.gradient.CL[bora::input_index::kAlpha] *
           kRadiansToDegrees;
}

Scalar BoraPolarSolver::GetDcdDalphaPerRadian() const {
    return impl_->result.gradient.CD[bora::input_index::kAlpha] *
           kRadiansToDegrees;
}

const Eigen::MatrixXd& BoraPolarSolver::GetDclDx() const {
    return impl_->dcl_dshape;
}

const Eigen::MatrixXd& BoraPolarSolver::GetDcdDx() const {
    return impl_->dcd_dshape;
}

Scalar BoraPolarSolver::GetAnalysisConfidence() const {
    return impl_->result.output.analysis_confidence;
}

std::vector<std::unique_ptr<AirfoilSolver>> CreateBoraStationSolvers(
    const std::vector<CouplingStation>& stations, const BoraOptions& options) {
    auto model = std::make_shared<const bora::NeuralFoil>(
        ToBoraModelSize(options.model_size),
        std::filesystem::path(options.data_dir));

    std::vector<std::unique_ptr<AirfoilSolver>> solvers;
    solvers.reserve(stations.size());
    for (const CouplingStation& station : stations) {
        auto solver = std::unique_ptr<BoraPolarSolver>(new BoraPolarSolver());
        solver->impl_->Initialize(model, station.coords, options);
        solvers.push_back(std::move(solver));
    }
    return solvers;
}

}  // namespace pegasus::viscous
