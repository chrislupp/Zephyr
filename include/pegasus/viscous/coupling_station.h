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

#ifndef PEGASUS_VISCOUS_COUPLING_STATION_H_
#define PEGASUS_VISCOUS_COUPLING_STATION_H_

#include "pegasus/core/types.h"
#include "pegasus/solver/vortex_lattice.h"

#include <vector>

#include <Eigen/Core>

namespace pegasus::viscous {

/// @brief A user-defined 2D viscous coupling station.
///
/// Stations are independent of the VLM strip grid: the user places
/// `N_stations` stations at chosen spanwise positions (eta values) with
/// their own airfoil coordinates and panel counts. The coupler interpolates
/// between strip data and station data using cubic splines along eta.
struct CouplingStation {
    /// Normalized spanwise position. The coupler defines eta = y / y_max,
    /// where y_max is the largest absolute spanwise coordinate in the mesh.
    /// For y-duplicated full wings, eta ranges over [−1, +1]; for half-
    /// wings it ranges over [0, 1].
    Scalar eta = 0.0;

    /// Airfoil node coordinates, shape (2, N). Passed directly to Scirocco.
    /// Order matches Scirocco's convention (clockwise from TE, sharp-TE,
    /// closed loop).
    Eigen::MatrixXd coords;

    /// Number of internal Scirocco panels (controls repaneling).
    int num_panels = 199;

    /// Per-station Reynolds number used at every Newton iteration.
    Scalar reynolds = 1.0e6;
};

/// @brief Strip↔station cubic-spline interpolation operators.
///
/// Built once at coupler construction time. Operates on the mesh strip
/// ordering (matching `mesh.GetStrip(i)`) and the user station ordering.
struct StationInterpolation {
    /// Strip eta values, one per VLM strip in mesh order.
    VecX eta_strips;

    /// Station eta values in user-supplied order (must be monotonically
    /// increasing).
    VecX eta_stations;

    /// Matrix mapping station values to strip values via cubic spline.
    /// Shape: (N_strips, N_stations).
    /// `strip_values = to_strip * station_values`
    MatX to_strip;

    /// Matrix mapping strip values to station values via cubic spline.
    /// Shape: (N_stations, N_strips).
    /// `station_values = to_station * strip_values`
    MatX to_station;
};

/// @brief Build the strip↔station spline operators for a given mesh and
/// station list.
///
/// Computes per-strip eta from `mesh.GetStrip(i).leading_edge_center.y()`
/// normalized by the maximum absolute spanwise extent in the mesh. The
/// `stations` vector must have at least 2 entries with strictly increasing
/// eta values, and the station eta range should bracket all strip etas
/// (otherwise the spline extrapolates and accuracy degrades near the tip).
///
/// Both operators are dense (cubic spline has global support through the
/// tridiagonal slope solve), so storage is O(N_strips · N_stations).
StationInterpolation BuildStationInterpolation(
    const VortexLattice& mesh, const std::vector<CouplingStation>& stations);

}  // namespace pegasus::viscous

#endif  // PEGASUS_VISCOUS_COUPLING_STATION_H_
