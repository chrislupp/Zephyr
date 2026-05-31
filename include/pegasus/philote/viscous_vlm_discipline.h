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
/// Inputs (matching VlmDiscipline's contract):
///   - One (3 * num_nodes) array per surface, named by the surface
///   - Reference quantities: reference_area / reference_chord /
///     reference_span / moment_ref_{x,y,z}
///   - Flight condition: alpha, beta, mach, dynamic_pressure,
///     roll_rate, pitch_rate, yaw_rate
///
/// Viscous options (via SetOptions):
///   - viscous_airfoil_naca (string, default "0012")
///   - viscous_reynolds (double, default 1e6)
///   - viscous_mach (double, default 0.0)
///   - viscous_panels (int, default 99)
///   - viscous_max_iter (int, default 10)
///   - viscous_tolerance (double, default 1e-4)
///
/// Outputs:
///   - CL, CD (viscous-corrected total coefficients)
///
/// Partials declared (v1):
///   - dCL/dalpha, dCD/dalpha (total derivatives via implicit-function
///     chain inside ViscousVlmCoupler::ComputeAlphaDerivatives())
///
/// Other partials (mesh nodes, ref quantities, other flight states) are
/// not yet implemented and return zero — to be added as the coupled
/// adjoint is extended in subsequent versions.
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
    viscous::CouplerOptions BuildCouplerOptions() const;

    std::vector<ViscousSurfaceConfig> surfaces_;

    // Coupler options (populated by SetOptions).
    std::string airfoil_naca_{"0012"};
    double reynolds_{1.0e6};
    double mach_override_{0.0};
    int scirocco_panels_{99};
    int max_iter_{10};
    double tolerance_{1.0e-4};
};

}  // namespace pegasus
