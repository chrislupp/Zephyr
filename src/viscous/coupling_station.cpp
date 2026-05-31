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

#include "pegasus/viscous/coupling_station.h"

#include "pegasus/utils/cubic_spline.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pegasus::viscous {

namespace {

/// Build a spline interpolation matrix: rows are evaluations, columns are
/// data points. Output(i, j) = spline(data_etas, e_j)(eval_etas[i]) where
/// e_j is the unit data vector. Uses natural cubic-spline boundary
/// conditions throughout.
MatX BuildSplineMatrix(const std::vector<Scalar>& data_etas,
                       const std::vector<Scalar>& eval_etas) {
    const Index n_data = static_cast<Index>(data_etas.size());
    const Index n_eval = static_cast<Index>(eval_etas.size());
    MatX M = MatX::Zero(n_eval, n_data);

    std::vector<Scalar> basis(static_cast<std::size_t>(n_data), 0.0);
    for (Index j = 0; j < n_data; ++j) {
        std::fill(basis.begin(), basis.end(), Scalar(0.0));
        basis[static_cast<std::size_t>(j)] = Scalar(1.0);
        CubicSpline spline;
        spline.Build(data_etas, basis);
        for (Index i = 0; i < n_eval; ++i) {
            M(i, j) = spline.Evaluate(eval_etas[static_cast<std::size_t>(i)]);
        }
    }
    return M;
}

}  // namespace

StationInterpolation BuildStationInterpolation(
    const VortexLattice& mesh, const std::vector<CouplingStation>& stations) {
    const Index n_strips = mesh.GetNumStrips();
    const Index n_stations = static_cast<Index>(stations.size());

    if (n_stations < 2) {
        throw std::invalid_argument(
            "BuildStationInterpolation: need at least 2 stations to build "
            "the cubic-spline operators.");
    }
    for (Index k = 1; k < n_stations; ++k) {
        if (stations[static_cast<std::size_t>(k)].eta <=
            stations[static_cast<std::size_t>(k - 1)].eta) {
            throw std::invalid_argument(
                "BuildStationInterpolation: station eta values must be "
                "strictly increasing.");
        }
    }

    // Strip etas: y / y_max where y_max is the largest absolute spanwise
    // coordinate over all strip leading-edge centers.
    Scalar y_max = 0.0;
    for (Index i = 0; i < n_strips; ++i) {
        Scalar y = mesh.GetStrip(i).leading_edge_center.y();
        if (std::abs(y) > y_max) {
            y_max = std::abs(y);
        }
    }
    if (y_max < Scalar(1e-30)) {
        throw std::invalid_argument(
            "BuildStationInterpolation: mesh has zero spanwise extent.");
    }

    StationInterpolation out;
    out.eta_strips.resize(n_strips);
    out.eta_stations.resize(n_stations);

    std::vector<Scalar> strip_eta_vec(static_cast<std::size_t>(n_strips));
    std::vector<Scalar> station_eta_vec(static_cast<std::size_t>(n_stations));
    for (Index i = 0; i < n_strips; ++i) {
        Scalar e = mesh.GetStrip(i).leading_edge_center.y() / y_max;
        out.eta_strips(i) = e;
        strip_eta_vec[static_cast<std::size_t>(i)] = e;
    }
    for (Index k = 0; k < n_stations; ++k) {
        Scalar e = stations[static_cast<std::size_t>(k)].eta;
        out.eta_stations(k) = e;
        station_eta_vec[static_cast<std::size_t>(k)] = e;
    }

    // station-data spline evaluated at strip etas: N_strips × N_stations.
    out.to_strip = BuildSplineMatrix(station_eta_vec, strip_eta_vec);

    // strip-data spline evaluated at station etas: N_stations × N_strips.
    // Cubic spline of strip-ordered data requires monotonic strip etas. The
    // mesh may not be naturally sorted that way (y-duplicate puts mirrored
    // strips after originals), so we sort indices by eta, build the spline
    // on the sorted data, then reorder the columns to match the original
    // mesh strip order.
    std::vector<Index> order(static_cast<std::size_t>(n_strips));
    for (Index i = 0; i < n_strips; ++i) {
        order[static_cast<std::size_t>(i)] = i;
    }
    std::sort(order.begin(), order.end(), [&](Index a, Index b) {
        return out.eta_strips(a) < out.eta_strips(b);
    });
    std::vector<Scalar> sorted_strip_eta(static_cast<std::size_t>(n_strips));
    for (Index i = 0; i < n_strips; ++i) {
        sorted_strip_eta[static_cast<std::size_t>(i)] =
            out.eta_strips(order[static_cast<std::size_t>(i)]);
    }
    // Guard against duplicate etas (which would break the tridiagonal solve).
    for (Index i = 1; i < n_strips; ++i) {
        if (sorted_strip_eta[static_cast<std::size_t>(i)] <=
            sorted_strip_eta[static_cast<std::size_t>(i - 1)]) {
            throw std::invalid_argument(
                "BuildStationInterpolation: mesh has duplicated strip eta "
                "values — spline cannot resolve.");
        }
    }

    MatX sorted_to_station =
        BuildSplineMatrix(sorted_strip_eta, station_eta_vec);
    out.to_station = MatX::Zero(n_stations, n_strips);
    for (Index i = 0; i < n_strips; ++i) {
        out.to_station.col(order[static_cast<std::size_t>(i)]) =
            sorted_to_station.col(i);
    }
    return out;
}

}  // namespace pegasus::viscous
