// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

#include "vlm/ad/gradient_calculator.h"
#include "vlm/aero/force_calculator.h"
#include "vlm/geometry/mesh_input.h"
#include "vlm/solver/vortex_lattice_solver.h"

#include <string>
#include <vector>

#include <explicit.h>

namespace vlm {

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
    ad::GradientCalculator grad_calc_;
};

}  // namespace vlm
