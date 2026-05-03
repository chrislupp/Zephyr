/*
    Project Pegasus
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

#include "pegasus/philote/vlm_discipline.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace pegasus {

// =============================================================================
// String parsing utilities
// =============================================================================

std::vector<std::string> VlmDiscipline::SplitString(const std::string& str,
                                                    char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(str);
    std::string token;
    while (std::getline(stream, token, delimiter)) {
        // Trim whitespace
        auto start = token.find_first_not_of(" \t");
        auto end = token.find_last_not_of(" \t");
        if (start != std::string::npos) {
            tokens.push_back(token.substr(start, end - start + 1));
        }
    }
    return tokens;
}

bool VlmDiscipline::ParseBool(const std::string& str) {
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower == "true" || lower == "1" || lower == "yes";
}

// =============================================================================
// Philote lifecycle methods
// =============================================================================

void VlmDiscipline::Initialize() {
    ExplicitDiscipline::Initialize();

    // Surface configuration
    AddOption("num_surfaces", "int");
    AddOption("surface_names", "string");
    AddOption("surface_nx", "string");
    AddOption("surface_ny", "string");
    AddOption("surface_has_wake", "string");
    AddOption("surface_y_duplicated", "string");
    AddOption("surface_y_duplicate_plane", "string");
    AddOption("surface_component_index", "string");
    AddOption("surface_affected_by_freestream", "string");
}

void VlmDiscipline::SetOptions(const google::protobuf::Struct& options_struct) {
    const auto& fields = options_struct.fields();

    // Parse number of surfaces
    int64_t num_surfaces = 1;
    auto it = fields.find("num_surfaces");
    if (it != fields.end()) {
        num_surfaces = static_cast<int64_t>(it->second.number_value());
    }

    // Parse surface names
    std::vector<std::string> names;
    it = fields.find("surface_names");
    if (it != fields.end()) {
        names = SplitString(it->second.string_value(), ',');
    }
    if (static_cast<int64_t>(names.size()) != num_surfaces) {
        throw std::runtime_error(
            "surface_names count does not match num_surfaces");
    }

    // Parse nx values
    std::vector<std::string> nx_strs;
    it = fields.find("surface_nx");
    if (it != fields.end()) {
        nx_strs = SplitString(it->second.string_value(), ',');
    }
    if (static_cast<int64_t>(nx_strs.size()) != num_surfaces) {
        throw std::runtime_error(
            "surface_nx count does not match num_surfaces");
    }

    // Parse ny values
    std::vector<std::string> ny_strs;
    it = fields.find("surface_ny");
    if (it != fields.end()) {
        ny_strs = SplitString(it->second.string_value(), ',');
    }
    if (static_cast<int64_t>(ny_strs.size()) != num_surfaces) {
        throw std::runtime_error(
            "surface_ny count does not match num_surfaces");
    }

    // Parse optional surface properties
    std::vector<std::string> has_wake_strs, y_dup_strs, y_dup_plane_strs,
        comp_idx_strs, aff_strs;

    it = fields.find("surface_has_wake");
    if (it != fields.end())
        has_wake_strs = SplitString(it->second.string_value(), ',');

    it = fields.find("surface_y_duplicated");
    if (it != fields.end())
        y_dup_strs = SplitString(it->second.string_value(), ',');

    it = fields.find("surface_y_duplicate_plane");
    if (it != fields.end())
        y_dup_plane_strs = SplitString(it->second.string_value(), ',');

    it = fields.find("surface_component_index");
    if (it != fields.end())
        comp_idx_strs = SplitString(it->second.string_value(), ',');

    it = fields.find("surface_affected_by_freestream");
    if (it != fields.end())
        aff_strs = SplitString(it->second.string_value(), ',');

    // Build surface configs
    surfaces_.clear();
    surfaces_.resize(num_surfaces);
    for (int64_t i = 0; i < num_surfaces; ++i) {
        surfaces_[i].name = names[i];
        surfaces_[i].nx = std::stoll(nx_strs[i]);
        surfaces_[i].ny = std::stoll(ny_strs[i]);

        if (i < static_cast<int64_t>(has_wake_strs.size()))
            surfaces_[i].has_wake = ParseBool(has_wake_strs[i]);
        if (i < static_cast<int64_t>(y_dup_strs.size()))
            surfaces_[i].y_duplicated = ParseBool(y_dup_strs[i]);
        if (i < static_cast<int64_t>(y_dup_plane_strs.size()))
            surfaces_[i].y_duplicate_plane = std::stod(y_dup_plane_strs[i]);
        if (i < static_cast<int64_t>(comp_idx_strs.size()))
            surfaces_[i].component_index = std::stoll(comp_idx_strs[i]);
        if (i < static_cast<int64_t>(aff_strs.size()))
            surfaces_[i].affected_by_freestream = ParseBool(aff_strs[i]);
    }

    // Call parent to trigger Configure()
    ExplicitDiscipline::SetOptions(options_struct);
}

void VlmDiscipline::Setup() {
    // Mesh node inputs (one per surface)
    for (const auto& s : surfaces_) {
        int64_t num_values = (s.nx + 1) * (s.ny + 1) * 3;
        AddInput(s.name, {num_values}, "m");
    }

    // Reference quantity inputs
    AddInput("reference_area", {1}, "m**2");
    AddInput("reference_chord", {1}, "m");
    AddInput("reference_span", {1}, "m");
    AddInput("moment_ref_x", {1}, "m");
    AddInput("moment_ref_y", {1}, "m");
    AddInput("moment_ref_z", {1}, "m");

    // Flight condition inputs
    AddInput("alpha", {1}, "rad");
    AddInput("beta", {1}, "rad");
    AddInput("mach", {1}, "");
    AddInput("dynamic_pressure", {1}, "Pa");
    AddInput("roll_rate", {1}, "");
    AddInput("pitch_rate", {1}, "");
    AddInput("yaw_rate", {1}, "");

    // Aerodynamic coefficient outputs
    AddOutput("CL", {1}, "");
    AddOutput("CD", {1}, "");
    AddOutput("CY", {1}, "");
    AddOutput("Cl", {1}, "");
    AddOutput("Cm", {1}, "");
    AddOutput("Cn", {1}, "");
}

void VlmDiscipline::SetupPartials() {
    const char* outputs[] = {"CL", "CD", "CY", "Cl", "Cm", "Cn"};

    // Partials w.r.t. mesh node coordinates
    for (const auto& s : surfaces_) {
        for (const char* out : outputs) {
            DeclarePartials(out, s.name);
        }
    }

    // Partials w.r.t. reference quantity inputs
    const char* ref_inputs[] = {"reference_area", "reference_chord",
                                "reference_span", "moment_ref_x",
                                "moment_ref_y",   "moment_ref_z"};
    for (const char* out : outputs) {
        for (const char* inp : ref_inputs) {
            DeclarePartials(out, inp);
        }
    }

    // Partials w.r.t. flight condition inputs
    const char* fc_inputs[] = {
        "alpha",     "beta",       "mach",    "dynamic_pressure",
        "roll_rate", "pitch_rate", "yaw_rate"};
    for (const char* out : outputs) {
        for (const char* inp : fc_inputs) {
            DeclarePartials(out, inp);
        }
    }
}

// =============================================================================
// Helper methods
// =============================================================================

MeshInput VlmDiscipline::BuildMeshInput(
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

FlightCondition VlmDiscipline::BuildFlightCondition(
    const philote::Variables& inputs) const {
    FlightCondition condition;
    condition.alpha = inputs.at("alpha")(0);
    condition.beta = inputs.at("beta")(0);
    condition.mach = inputs.at("mach")(0);
    condition.dynamic_pressure = inputs.at("dynamic_pressure")(0);
    condition.roll_rate = inputs.at("roll_rate")(0);
    condition.pitch_rate = inputs.at("pitch_rate")(0);
    condition.yaw_rate = inputs.at("yaw_rate")(0);
    return condition;
}

// =============================================================================
// Compute
// =============================================================================

void VlmDiscipline::Compute(const philote::Variables& inputs,
                            philote::Variables& outputs) {
    MeshInput mesh_input = BuildMeshInput(inputs);
    FlightCondition condition = BuildFlightCondition(inputs);

    solver_.SetMesh(mesh_input);
    solver_.RebuildAIC();
    SolutionState state = solver_.Solve(condition);

    ForceCalculator calc;
    AerodynamicResults results = calc.ComputeForces(
        mesh_input.reference, solver_.GetMesh(), state, condition);

    outputs.at("CL")(0) = results.CL;
    outputs.at("CD")(0) = results.CD;
    outputs.at("CY")(0) = results.CY;
    outputs.at("Cl")(0) = results.Cl;
    outputs.at("Cm")(0) = results.Cm;
    outputs.at("Cn")(0) = results.Cn;
}

// =============================================================================
// ComputePartials
// =============================================================================

void VlmDiscipline::ComputePartials(const philote::Variables& inputs,
                                    philote::Partials& partials) {
    MeshInput mesh_input = BuildMeshInput(inputs);
    FlightCondition condition = BuildFlightCondition(inputs);

    // We need the coefficient values for reference quantity partials
    grad_calc_.SetMesh(mesh_input);

    // Solve forward problem (also caches for gradient computation)
    solver_.SetMesh(mesh_input);
    solver_.RebuildAIC();
    SolutionState state = solver_.Solve(condition);
    ForceCalculator calc;
    AerodynamicResults results = calc.ComputeForces(
        mesh_input.reference, solver_.GetMesh(), state, condition);

    // Mesh node gradients
    auto node_grads = grad_calc_.ComputeMeshNodeGradients(condition);

    for (size_t s = 0; s < surfaces_.size(); ++s) {
        const auto& sg = node_grads[s];
        const std::string& name = surfaces_[s].name;
        int64_t n = (surfaces_[s].nx + 1) * (surfaces_[s].ny + 1) * 3;

        for (int64_t i = 0; i < n; ++i) {
            partials[{"CL", name}](i) = sg.dCL.d_nodes[i];
            partials[{"CD", name}](i) = sg.dCD.d_nodes[i];
            partials[{"CY", name}](i) = sg.dCY.d_nodes[i];
            partials[{"Cl", name}](i) = sg.dCl.d_nodes[i];
            partials[{"Cm", name}](i) = sg.dCm.d_nodes[i];
            partials[{"Cn", name}](i) = sg.dCn.d_nodes[i];
        }
    }

    // Flight condition gradients
    gradients::FullGradients fc_grads =
        grad_calc_.ComputeFlightConditionDerivatives(condition);

    // CL partials
    partials[{"CL", "alpha"}](0) = fc_grads.alpha_derivs.dCL;
    partials[{"CL", "beta"}](0) = fc_grads.beta_derivs.dCL;
    partials[{"CL", "mach"}](0) = fc_grads.mach_derivs.dCL;
    partials[{"CL", "roll_rate"}](0) = fc_grads.p_derivs.dCL;
    partials[{"CL", "pitch_rate"}](0) = fc_grads.q_derivs.dCL;
    partials[{"CL", "yaw_rate"}](0) = fc_grads.r_derivs.dCL;
    partials[{"CL", "dynamic_pressure"}](0) = 0.0;

    // CD partials
    partials[{"CD", "alpha"}](0) = fc_grads.alpha_derivs.dCD;
    partials[{"CD", "beta"}](0) = fc_grads.beta_derivs.dCD;
    partials[{"CD", "mach"}](0) = fc_grads.mach_derivs.dCD;
    partials[{"CD", "roll_rate"}](0) = fc_grads.p_derivs.dCD;
    partials[{"CD", "pitch_rate"}](0) = fc_grads.q_derivs.dCD;
    partials[{"CD", "yaw_rate"}](0) = fc_grads.r_derivs.dCD;
    partials[{"CD", "dynamic_pressure"}](0) = 0.0;

    // CY partials
    partials[{"CY", "alpha"}](0) = fc_grads.alpha_derivs.dCY;
    partials[{"CY", "beta"}](0) = fc_grads.beta_derivs.dCY;
    partials[{"CY", "mach"}](0) = fc_grads.mach_derivs.dCY;
    partials[{"CY", "roll_rate"}](0) = fc_grads.p_derivs.dCY;
    partials[{"CY", "pitch_rate"}](0) = fc_grads.q_derivs.dCY;
    partials[{"CY", "yaw_rate"}](0) = fc_grads.r_derivs.dCY;
    partials[{"CY", "dynamic_pressure"}](0) = 0.0;

    // Cl partials
    partials[{"Cl", "alpha"}](0) = fc_grads.alpha_derivs.dCl;
    partials[{"Cl", "beta"}](0) = fc_grads.beta_derivs.dCl;
    partials[{"Cl", "mach"}](0) = fc_grads.mach_derivs.dCl;
    partials[{"Cl", "roll_rate"}](0) = fc_grads.p_derivs.dCl;
    partials[{"Cl", "pitch_rate"}](0) = fc_grads.q_derivs.dCl;
    partials[{"Cl", "yaw_rate"}](0) = fc_grads.r_derivs.dCl;
    partials[{"Cl", "dynamic_pressure"}](0) = 0.0;

    // Cm partials
    partials[{"Cm", "alpha"}](0) = fc_grads.alpha_derivs.dCm;
    partials[{"Cm", "beta"}](0) = fc_grads.beta_derivs.dCm;
    partials[{"Cm", "mach"}](0) = fc_grads.mach_derivs.dCm;
    partials[{"Cm", "roll_rate"}](0) = fc_grads.p_derivs.dCm;
    partials[{"Cm", "pitch_rate"}](0) = fc_grads.q_derivs.dCm;
    partials[{"Cm", "yaw_rate"}](0) = fc_grads.r_derivs.dCm;
    partials[{"Cm", "dynamic_pressure"}](0) = 0.0;

    // Cn partials
    partials[{"Cn", "alpha"}](0) = fc_grads.alpha_derivs.dCn;
    partials[{"Cn", "beta"}](0) = fc_grads.beta_derivs.dCn;
    partials[{"Cn", "mach"}](0) = fc_grads.mach_derivs.dCn;
    partials[{"Cn", "roll_rate"}](0) = fc_grads.p_derivs.dCn;
    partials[{"Cn", "pitch_rate"}](0) = fc_grads.q_derivs.dCn;
    partials[{"Cn", "yaw_rate"}](0) = fc_grads.r_derivs.dCn;
    partials[{"Cn", "dynamic_pressure"}](0) = 0.0;

    // Reference quantity partials
    // All coefficients are proportional to 1/S_ref:
    //   d(coeff)/d(S_ref) = -coeff / S_ref
    double s_ref = mesh_input.reference.reference_area;
    double c_ref = mesh_input.reference.reference_chord;

    // d/d(reference_area): all coefficients scale as 1/S_ref
    partials[{"CL", "reference_area"}](0) = -results.CL / s_ref;
    partials[{"CD", "reference_area"}](0) = -results.CD / s_ref;
    partials[{"CY", "reference_area"}](0) = -results.CY / s_ref;
    partials[{"Cl", "reference_area"}](0) = -results.Cl / s_ref;
    partials[{"Cm", "reference_area"}](0) = -results.Cm / s_ref;
    partials[{"Cn", "reference_area"}](0) = -results.Cn / s_ref;

    // d/d(reference_chord): moment coefficients also scale as 1/c_ref
    partials[{"CL", "reference_chord"}](0) = 0.0;
    partials[{"CD", "reference_chord"}](0) = 0.0;
    partials[{"CY", "reference_chord"}](0) = 0.0;
    partials[{"Cl", "reference_chord"}](0) = -results.Cl / c_ref;
    partials[{"Cm", "reference_chord"}](0) = -results.Cm / c_ref;
    partials[{"Cn", "reference_chord"}](0) = -results.Cn / c_ref;

    // d/d(reference_span): force/moment coefficients don't depend on b_ref
    partials[{"CL", "reference_span"}](0) = 0.0;
    partials[{"CD", "reference_span"}](0) = 0.0;
    partials[{"CY", "reference_span"}](0) = 0.0;
    partials[{"Cl", "reference_span"}](0) = 0.0;
    partials[{"Cm", "reference_span"}](0) = 0.0;
    partials[{"Cn", "reference_span"}](0) = 0.0;

    // d/d(moment_ref_x/y/z): zero for now (moment ref point affects moments
    // but computing the exact partial requires the total force vector)
    const char* moment_ref_inputs[] = {"moment_ref_x", "moment_ref_y",
                                       "moment_ref_z"};
    const char* all_outputs[] = {"CL", "CD", "CY", "Cl", "Cm", "Cn"};
    for (const char* inp : moment_ref_inputs) {
        for (const char* out : all_outputs) {
            partials[{out, inp}](0) = 0.0;
        }
    }
}

}  // namespace pegasus
