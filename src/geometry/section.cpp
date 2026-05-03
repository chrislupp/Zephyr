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

#include "pegasus/geometry/section.h"

#include <utility>

namespace pegasus {

namespace {
// Static flat plate airfoil for sections without custom airfoil
Airfoil g_flat_plate = Airfoil::CreateFlatPlate();
}  // namespace

Section::Section(Vec3 leading_edge_position, Scalar chord, Scalar twist_angle,
                 std::shared_ptr<Airfoil> airfoil)
    : leading_edge_position_(std::move(leading_edge_position)),
      chord_(chord),
      twist_angle_deg_(twist_angle),
      airfoil_(std::move(std::move(airfoil))) {}

const Airfoil& Section::GetAirfoil() const {
    if (airfoil_) {
        return *airfoil_;
    }
    return g_flat_plate;
}

Airfoil& Section::GetAirfoilMutable() {
    if (!airfoil_) {
        airfoil_ = std::make_shared<Airfoil>(Airfoil::CreateFlatPlate());
    }
    return *airfoil_;
}

}  // namespace pegasus
