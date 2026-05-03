// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/geometry/section.h"

namespace vlm {

namespace {
// Static flat plate airfoil for sections without custom airfoil
Airfoil g_flat_plate = Airfoil::CreateFlatPlate();
}  // namespace

Section::Section(const Vec3& leading_edge_position, Scalar chord,
                 Scalar twist_angle, std::shared_ptr<Airfoil> airfoil)
    : leading_edge_position_(leading_edge_position),
      chord_(chord),
      twist_angle_deg_(twist_angle),
      airfoil_(airfoil) {}

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

}  // namespace vlm
