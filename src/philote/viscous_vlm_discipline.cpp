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

#include "pegasus/philote/viscous_vlm_discipline.h"

#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/core/types.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pegasus {

namespace {

std::vector<std::string> SplitCsv(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream stream(str);
    std::string token;
    while (std::getline(stream, token, ',')) {
        auto start = token.find_first_not_of(" \t");
        auto end = token.find_last_not_of(" \t");
        if (start != std::string::npos) {
            tokens.push_back(token.substr(start, end - start + 1));
        }
    }
    return tokens;
}

bool ParseBoolStr(const std::string& s) {
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower == "true" || lower == "1" || lower == "yes";
}

std::string StationCoordsInputName(int k) {
    std::ostringstream oss;
    oss << "station_" << k << "_coords";
    return oss.str();
}

}  // namespace

// =============================================================================
// Philote lifecycle
// =============================================================================

void ViscousVlmDiscipline::Initialize() {
    ExplicitDiscipline::Initialize();

    // Surface configuration.
    AddOption("num_surfaces", "int");
    AddOption("surface_names", "string");
    AddOption("surface_nx", "string");
    AddOption("surface_ny", "string");
    AddOption("surface_has_wake", "string");
    AddOption("surface_y_duplicated", "string");
    AddOption("surface_y_duplicate_plane", "string");
    AddOption("surface_component_index", "string");
    AddOption("surface_affected_by_freestream", "string");

    // Per-station configuration.
    AddOption("viscous_num_stations", "int");
    AddOption("viscous_station_eta", "string");
    AddOption("viscous_station_panels", "string");
    AddOption("viscous_station_reynolds", "string");

    // Global viscous-coupling scalars.
    AddOption("viscous_mach", "double");
    AddOption("viscous_max_iter", "int");
    AddOption("viscous_tolerance", "double");
    AddOption("viscous_max_step_deg", "double");
}

void ViscousVlmDiscipline::SetOptions(
    const google::protobuf::Struct& options_struct) {
    const auto& fields = options_struct.fields();

    auto read_csv = [&](const std::string& key) {
        auto f = fields.find(key);
        return f != fields.end() ? SplitCsv(f->second.string_value())
                                 : std::vector<std::string>{};
    };

    // -------- surfaces (same as VlmDiscipline) --------
    int64_t num_surfaces = 1;
    auto it = fields.find("num_surfaces");
    if (it != fields.end()) {
        num_surfaces = static_cast<int64_t>(it->second.number_value());
    }
    auto names = read_csv("surface_names");
    auto nx_strs = read_csv("surface_nx");
    auto ny_strs = read_csv("surface_ny");
    auto has_wake_strs = read_csv("surface_has_wake");
    auto y_dup_strs = read_csv("surface_y_duplicated");
    auto y_dup_plane_strs = read_csv("surface_y_duplicate_plane");
    auto comp_idx_strs = read_csv("surface_component_index");
    auto aff_strs = read_csv("surface_affected_by_freestream");
    if (static_cast<int64_t>(names.size()) != num_surfaces ||
        static_cast<int64_t>(nx_strs.size()) != num_surfaces ||
        static_cast<int64_t>(ny_strs.size()) != num_surfaces) {
        throw std::runtime_error(
            "ViscousVlmDiscipline: surface_names / surface_nx / surface_ny "
            "count must match num_surfaces");
    }
    surfaces_.clear();
    surfaces_.resize(num_surfaces);
    for (int64_t i = 0; i < num_surfaces; ++i) {
        surfaces_[i].name = names[i];
        surfaces_[i].nx = std::stoll(nx_strs[i]);
        surfaces_[i].ny = std::stoll(ny_strs[i]);
        if (i < static_cast<int64_t>(has_wake_strs.size()))
            surfaces_[i].has_wake = ParseBoolStr(has_wake_strs[i]);
        if (i < static_cast<int64_t>(y_dup_strs.size()))
            surfaces_[i].y_duplicated = ParseBoolStr(y_dup_strs[i]);
        if (i < static_cast<int64_t>(y_dup_plane_strs.size()))
            surfaces_[i].y_duplicate_plane = std::stod(y_dup_plane_strs[i]);
        if (i < static_cast<int64_t>(comp_idx_strs.size()))
            surfaces_[i].component_index = std::stoll(comp_idx_strs[i]);
        if (i < static_cast<int64_t>(aff_strs.size()))
            surfaces_[i].affected_by_freestream = ParseBoolStr(aff_strs[i]);
    }

    // -------- stations --------
    int64_t num_stations = 0;
    it = fields.find("viscous_num_stations");
    if (it != fields.end()) {
        num_stations = static_cast<int64_t>(it->second.number_value());
    }
    if (num_stations < 2) {
        throw std::runtime_error(
            "ViscousVlmDiscipline: viscous_num_stations must be >= 2");
    }
    auto eta_strs = read_csv("viscous_station_eta");
    auto panel_strs = read_csv("viscous_station_panels");
    auto re_strs = read_csv("viscous_station_reynolds");
    if (static_cast<int64_t>(eta_strs.size()) != num_stations) {
        throw std::runtime_error(
            "ViscousVlmDiscipline: viscous_station_eta count must match "
            "viscous_num_stations");
    }
    station_eta_.resize(num_stations);
    station_panels_.assign(num_stations, 99);
    station_reynolds_.assign(num_stations, 1.0e6);
    for (int64_t i = 0; i < num_stations; ++i) {
        station_eta_[i] = std::stod(eta_strs[i]);
        if (i < static_cast<int64_t>(panel_strs.size())) {
            station_panels_[i] = static_cast<int>(std::stol(panel_strs[i]));
        }
        if (i < static_cast<int64_t>(re_strs.size())) {
            station_reynolds_[i] = std::stod(re_strs[i]);
        }
    }

    // -------- scalars --------
    it = fields.find("viscous_mach");
    if (it != fields.end())
        mach_override_ = it->second.number_value();
    it = fields.find("viscous_max_iter");
    if (it != fields.end())
        max_iter_ = static_cast<int>(it->second.number_value());
    it = fields.find("viscous_tolerance");
    if (it != fields.end())
        tolerance_ = it->second.number_value();
    it = fields.find("viscous_max_step_deg");
    if (it != fields.end())
        max_step_deg_ = it->second.number_value();

    ExplicitDiscipline::SetOptions(options_struct);
}

void ViscousVlmDiscipline::Setup() {
    // Mesh node inputs (one per surface).
    for (const auto& s : surfaces_) {
        int64_t num_values = (s.nx + 1) * (s.ny + 1) * 3;
        AddInput(s.name, {num_values}, "m");
    }
    AddInput("reference_area", {1}, "m**2");
    AddInput("reference_chord", {1}, "m");
    AddInput("reference_span", {1}, "m");
    AddInput("moment_ref_x", {1}, "m");
    AddInput("moment_ref_y", {1}, "m");
    AddInput("moment_ref_z", {1}, "m");
    AddInput("alpha", {1}, "rad");
    AddInput("beta", {1}, "rad");
    AddInput("mach", {1}, "");
    AddInput("dynamic_pressure", {1}, "Pa");
    AddInput("roll_rate", {1}, "");
    AddInput("pitch_rate", {1}, "");
    AddInput("yaw_rate", {1}, "");

    // Per-station coords: 2 · N_k coords stored as [x0,z0,x1,z1,...].
    for (std::size_t k = 0; k < station_eta_.size(); ++k) {
        int n = station_panels_[k];
        AddInput(StationCoordsInputName(static_cast<int>(k)),
                 {2L * static_cast<int64_t>(n)}, "m");
    }

    AddOutput("CL", {1}, "");
    AddOutput("CD", {1}, "");
}

void ViscousVlmDiscipline::SetupPartials() {
    DeclarePartials("CL", "alpha");
    DeclarePartials("CD", "alpha");
    for (std::size_t k = 0; k < station_eta_.size(); ++k) {
        std::string name = StationCoordsInputName(static_cast<int>(k));
        DeclarePartials("CL", name);
        DeclarePartials("CD", name);
    }
}

// =============================================================================
// Helpers
// =============================================================================

MeshInput ViscousVlmDiscipline::BuildMeshInput(
    const philote::Variables& inputs) const {
    MeshInput mesh_input;
    mesh_input.reference.reference_area = inputs.at("reference_area")(0);
    mesh_input.reference.reference_chord = inputs.at("reference_chord")(0);
    mesh_input.reference.reference_span = inputs.at("reference_span")(0);
    mesh_input.reference.moment_reference_point[0] =
        inputs.at("moment_ref_x")(0);
    mesh_input.reference.moment_reference_point[1] =
        inputs.at("moment_ref_y")(0);
    mesh_input.reference.moment_reference_point[2] =
        inputs.at("moment_ref_z")(0);
    for (const auto& s : surfaces_) {
        SurfaceMesh surface;
        surface.name = s.name;
        surface.nx = s.nx;
        surface.ny = s.ny;
        surface.component_index = s.component_index;
        surface.affected_by_freestream = s.affected_by_freestream;
        surface.has_wake = s.has_wake;
        surface.y_duplicated = s.y_duplicated;
        surface.y_duplicate_plane = s.y_duplicate_plane;
        int64_t num_nodes = (s.nx + 1) * (s.ny + 1);
        surface.nodes.resize(num_nodes, 3);
        const auto& var = inputs.at(s.name);
        for (int64_t n = 0; n < num_nodes; ++n) {
            surface.nodes(n, 0) = var(n * 3 + 0);
            surface.nodes(n, 1) = var(n * 3 + 1);
            surface.nodes(n, 2) = var(n * 3 + 2);
        }
        mesh_input.surfaces.push_back(std::move(surface));
    }
    return mesh_input;
}

FlightCondition ViscousVlmDiscipline::BuildFlightCondition(
    const philote::Variables& inputs) const {
    FlightCondition c;
    c.alpha = inputs.at("alpha")(0);
    c.beta = inputs.at("beta")(0);
    c.mach = inputs.at("mach")(0);
    c.dynamic_pressure = inputs.at("dynamic_pressure")(0);
    c.roll_rate = inputs.at("roll_rate")(0);
    c.pitch_rate = inputs.at("pitch_rate")(0);
    c.yaw_rate = inputs.at("yaw_rate")(0);
    return c;
}

std::vector<viscous::CouplingStation> ViscousVlmDiscipline::BuildStations(
    const philote::Variables& inputs) const {
    std::vector<viscous::CouplingStation> stations;
    stations.reserve(station_eta_.size());
    for (std::size_t k = 0; k < station_eta_.size(); ++k) {
        viscous::CouplingStation s;
        s.eta = station_eta_[k];
        s.num_panels = station_panels_[k];
        s.reynolds = station_reynolds_[k];
        const auto& var =
            inputs.at(StationCoordsInputName(static_cast<int>(k)));
        const int N = station_panels_[k];
        s.coords.resize(2, N);
        for (int j = 0; j < N; ++j) {
            s.coords(0, j) = var(2 * j + 0);
            s.coords(1, j) = var(2 * j + 1);
        }
        stations.push_back(std::move(s));
    }
    return stations;
}

viscous::CouplerOptions ViscousVlmDiscipline::BuildCouplerOptions(
    const philote::Variables& inputs) const {
    viscous::CouplerOptions opt;
    opt.stations = BuildStations(inputs);
    opt.mach = mach_override_;
    opt.max_iterations = max_iter_;
    opt.tolerance = tolerance_;
    opt.max_step_rad = max_step_deg_ * kDegreesToRadians;
    return opt;
}

// =============================================================================
// Compute / ComputePartials
// =============================================================================

void ViscousVlmDiscipline::Compute(const philote::Variables& inputs,
                                   philote::Variables& outputs) {
    MeshInput mesh_input = BuildMeshInput(inputs);
    FlightCondition condition = BuildFlightCondition(inputs);
    VortexLatticeSolver solver;
    solver.SetMesh(mesh_input);
    viscous::ViscousVlmCoupler coupler(solver, condition,
                                       BuildCouplerOptions(inputs));
    viscous::CouplerResult res = coupler.Solve();
    outputs.at("CL")(0) = res.results.CL;
    outputs.at("CD")(0) = res.results.CD;
}

void ViscousVlmDiscipline::ComputePartials(const philote::Variables& inputs,
                                           philote::Partials& partials) {
    MeshInput mesh_input = BuildMeshInput(inputs);
    FlightCondition condition = BuildFlightCondition(inputs);
    VortexLatticeSolver solver;
    solver.SetMesh(mesh_input);
    viscous::ViscousVlmCoupler coupler(solver, condition,
                                       BuildCouplerOptions(inputs));
    viscous::CouplerResult res = coupler.Solve();
    if (!res.converged) {
        throw std::runtime_error(
            "ViscousVlmDiscipline::ComputePartials: viscous coupling did not "
            "converge");
    }

    viscous::CoupledGradients g_alpha = coupler.ComputeAlphaDerivatives();
    partials.at({"CL", "alpha"})(0) = g_alpha.dCL_dalpha;
    partials.at({"CD", "alpha"})(0) = g_alpha.dCD_dalpha;

    viscous::CoupledShapeGradients g_shape = coupler.ComputeShapeDerivatives();
    for (std::size_t k = 0; k < g_shape.dCL_dcoords.size(); ++k) {
        const Eigen::MatrixXd& dCL_dx = g_shape.dCL_dcoords[k];
        const Eigen::MatrixXd& dCD_dx = g_shape.dCD_dcoords[k];
        std::string name = StationCoordsInputName(static_cast<int>(k));
        auto& cl_p = partials.at({"CL", name});
        auto& cd_p = partials.at({"CD", name});
        const int N = static_cast<int>(dCL_dx.cols());
        // Coords are flattened as [x0,z0,x1,z1,…]; partials match.
        for (int j = 0; j < N; ++j) {
            cl_p(2 * j + 0) = dCL_dx(0, j);
            cl_p(2 * j + 1) = dCL_dx(1, j);
            cd_p(2 * j + 0) = dCD_dx(0, j);
            cd_p(2 * j + 1) = dCD_dx(1, j);
        }
    }
}

}  // namespace pegasus
