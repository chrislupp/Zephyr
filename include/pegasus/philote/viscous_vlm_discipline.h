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

#pragma once

#include "pegasus/geometry/mesh_input.h"
#include "pegasus/solver/vortex_lattice_solver.h"
#include "pegasus/viscous/coupler.h"

#include <string>
#include <vector>

#include <explicit.h>

namespace pegasus {

/// @brief Per-surface mesh config used by ViscousVlmDiscipline.
struct ViscousSurfaceConfig {
    std::string name;
    int64_t nx{0};
    int64_t ny{0};
    int64_t component_index{0};
    bool affected_by_freestream{true};
    bool has_wake{true};
    bool y_duplicated{false};
    double y_duplicate_plane{0.0};
};

/// @brief Philote discipline exposing the viscous-corrected VLM ↔ 2D
/// coupled aerodynamic outputs as a single monolithic MDO component.
///
/// Mesh / flight-condition inputs (matching VlmDiscipline):
///   - One (3 * num_nodes) array per surface, named by the surface
///   - Reference quantities: reference_area / reference_chord /
///     reference_span / moment_ref_{x,y,z}
///   - Flight condition: alpha, beta, mach, dynamic_pressure,
///     roll_rate, pitch_rate, yaw_rate
///
/// Per-station coupling inputs (one set per station, indexed by k):
///   - `station_k_coords` : flattened 2·N_k vector of airfoil node
///                          coordinates in Scirocco's order, stored
///                          row-major as [x0, z0, x1, z1, …].
///
/// Options (via SetOptions):
///   - Surface options (num_surfaces, surface_names, surface_nx, …)
///   - viscous_num_stations (int)
///   - viscous_station_eta  (CSV of doubles, length = num_stations)
///   - viscous_station_panels (CSV of ints; default 99 each if empty)
///   - viscous_station_reynolds (CSV of doubles; default 1e6 each if empty)
///   - viscous_mach (double, default 0.0)
///   - viscous_max_iter (int, default 30)
///   - viscous_tolerance (double, default 1e-6)
///   - viscous_max_step_deg (double, default 1.0)
///
/// Outputs:
///   - CL, CD (viscous-corrected total coefficients)
///
/// Partials declared:
///   - dCL/dalpha, dCD/dalpha (total derivatives via
///     ComputeAlphaDerivatives)
///   - dCL/dstation_k_coords, dCD/dstation_k_coords for each station k
///     (total derivatives via ComputeShapeDerivatives)
///
/// Other partials (mesh nodes, ref quantities, β, Mach, rotation rates)
/// are not yet implemented and return zero.
class ViscousVlmDiscipline : public philote::ExplicitDiscipline {
public:
    ViscousVlmDiscipline() = default;
    ~ViscousVlmDiscipline() override = default;

private:
    void Initialize() override;
    void SetOptions(const google::protobuf::Struct& options_struct) override;
    void Setup() override;
    void SetupPartials() override;

    void Compute(const philote::Variables& inputs,
                 philote::Variables& outputs) override;

    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials& partials) override;

    MeshInput BuildMeshInput(const philote::Variables& inputs) const;
    FlightCondition BuildFlightCondition(
        const philote::Variables& inputs) const;
    /// Build the station list by combining the per-station options (eta,
    /// panels, reynolds) with the per-station `station_k_coords` input.
    std::vector<viscous::CouplingStation> BuildStations(
        const philote::Variables& inputs) const;
    viscous::CouplerOptions BuildCouplerOptions(
        const philote::Variables& inputs) const;

    std::vector<ViscousSurfaceConfig> surfaces_;

    // Per-station configuration (populated by SetOptions).
    std::vector<double> station_eta_;
    std::vector<int> station_panels_;
    std::vector<double> station_reynolds_;

    // Global viscous-coupling scalars.
    double mach_override_{0.0};
    int max_iter_{30};
    double tolerance_{1.0e-6};
    double max_step_deg_{1.0};
};

}  // namespace pegasus
