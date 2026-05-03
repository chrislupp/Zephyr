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

#include "pegasus/geometry/surface.h"

#include <cmath>
#include <utility>

namespace pegasus {

Surface::Surface(std::string name) : name_(std::move(name)) {}

bool Surface::HasViscousDrag() const {
    // Check if any polar coefficient is non-zero
    for (const auto& coeff : cl_cd_polar_) {
        if (std::abs(coeff) > kEpsilon) {
            return true;
        }
    }
    return false;
}

Scalar Surface::ComputePlanformArea() const {
    if (sections_.size() < 2) {
        return 0.0;
    }

    Scalar area = 0.0;
    for (size_t i = 0; i < sections_.size() - 1; ++i) {
        const Section& sec1 = sections_[i];
        const Section& sec2 = sections_[i + 1];

        // Trapezoidal area between sections
        Scalar avg_chord = 0.5 * (sec1.GetChord() + sec2.GetChord());
        Vec3 delta =
            sec2.GetLeadingEdgePosition() - sec1.GetLeadingEdgePosition();
        Scalar span = std::sqrt(delta.y() * delta.y() + delta.z() * delta.z());

        area += avg_chord * span;
    }

    // Double if Y-duplicated
    if (IsYDuplicated()) {
        area *= 2.0;
    }

    return area;
}

Scalar Surface::ComputeMeanAerodynamicChord() const {
    if (sections_.size() < 2) {
        return sections_.empty() ? 0.0 : sections_[0].GetChord();
    }

    // MAC = integral(c^2 dy) / integral(c dy)
    // Using trapezoidal integration
    Scalar num = 0.0;
    Scalar den = 0.0;

    for (size_t i = 0; i < sections_.size() - 1; ++i) {
        const Section& sec1 = sections_[i];
        const Section& sec2 = sections_[i + 1];

        Scalar c1 = sec1.GetChord();
        Scalar c2 = sec2.GetChord();
        Vec3 delta =
            sec2.GetLeadingEdgePosition() - sec1.GetLeadingEdgePosition();
        Scalar dy = std::sqrt(delta.y() * delta.y() + delta.z() * delta.z());

        // Trapezoidal integration
        num += 0.5 * (c1 * c1 + c2 * c2) * dy;
        den += 0.5 * (c1 + c2) * dy;
    }

    return (den > kEpsilon) ? (num / den) : 0.0;
}

Scalar Surface::ComputeSpan() const {
    if (sections_.empty()) {
        return 0.0;
    }

    Scalar y_min = sections_[0].GetLeadingEdgePosition().y();
    Scalar y_max = y_min;

    for (const auto& section : sections_) {
        Scalar y = section.GetLeadingEdgePosition().y();
        y_min = std::min(y_min, y);
        y_max = std::max(y_max, y);
    }

    Scalar span = y_max - y_min;

    // Double if Y-duplicated
    if (IsYDuplicated()) {
        span *= 2.0;
    }

    return span;
}

}  // namespace pegasus
