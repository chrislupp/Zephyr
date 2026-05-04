/*
    Project Pegasus
    Version 0.1.0

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

#include <cstdint>
#include <string>
#include <vector>

#include <explicit.h>

namespace pegasus {

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
    bool y_duplicated_ = true;

    // Precomputed chordwise fractions [nx+1]
    std::vector<double> chord_fracs_;
};

}  // namespace pegasus
