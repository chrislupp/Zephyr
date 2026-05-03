// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <explicit.h>

namespace vlm {

/// Philote discipline that generates a structured surface mesh from
/// per-spanwise-station geometric parameters (chord, twist, leading-edge
/// position). Follows the OpenAeroStruct convention where design variables
/// are defined at each spanwise node directly.
///
/// Options define the panel counts and spacing type. Inputs are arrays of
/// size ny+1 (one per spanwise station). The output is a flat mesh node
/// array compatible with VlmDiscipline's surface input.
class MeshDiscipline : public philote::ExplicitDiscipline {
public:
    MeshDiscipline() = default;
    ~MeshDiscipline() = default;

private:
    void Initialize() override;
    void SetOptions(const google::protobuf::Struct& options_struct) override;
    void Setup() override;
    void SetupPartials() override;

    void Compute(const philote::Variables& inputs,
                 philote::Variables& outputs) override;

    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials& partials) override;

    void PrecomputeChordFracs();

    // Options
    int64_t nx_ = 1;
    int64_t ny_ = 1;
    bool cosine_chordwise_spacing_ = false;
    bool cosine_spanwise_spacing_ = false;

    // Precomputed chordwise fractions [nx+1]
    std::vector<double> chord_fracs_;
};

}  // namespace vlm
