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

#pragma once

#include "pegasus/gradients/gradient_calculator.h"
#include "pegasus/aero/force_calculator.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <string>
#include <vector>

#include <explicit.h>

namespace pegasus {

struct SurfaceConfig {
    std::string name;
    int64_t nx;
    int64_t ny;
    int64_t component_index = 0;
    bool affected_by_freestream = true;
    bool has_wake = true;
    bool y_duplicated = false;
    double y_duplicate_plane = 0.0;
};

class VlmDiscipline : public philote::ExplicitDiscipline {
public:
    VlmDiscipline() = default;
    ~VlmDiscipline() = default;

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

    static std::vector<std::string> SplitString(const std::string& str,
                                                char delimiter);
    static bool ParseBool(const std::string& str);

    std::vector<SurfaceConfig> surfaces_;
    VortexLatticeSolver solver_;
    gradients::GradientCalculator grad_calc_;
};

}  // namespace pegasus
