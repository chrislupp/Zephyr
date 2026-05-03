/*
    Project VLM
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

#include "vlm/geometry/aircraft.h"

#include <algorithm>
#include <set>

namespace vlm {

Aircraft::Aircraft(const std::string& title) : title_(title) {}

const Surface* Aircraft::FindSurfaceByName(const std::string& name) const {
    for (const auto& surface : surfaces_) {
        if (surface.GetName() == name) {
            return &surface;
        }
    }
    return nullptr;
}

std::vector<std::string> Aircraft::GetControlSurfaceNames() const {
    std::set<std::string> names;
    for (const auto& surface : surfaces_) {
        for (const auto& section : surface.GetSections()) {
            for (const auto& control : section.GetControlSurfaces()) {
                names.insert(control.name);
            }
        }
    }
    return std::vector<std::string>(names.begin(), names.end());
}

size_t Aircraft::GetNumControlSurfaces() const {
    return GetControlSurfaceNames().size();
}

std::vector<std::string> Aircraft::GetDesignVariableNames() const {
    std::set<std::string> names;
    for (const auto& surface : surfaces_) {
        for (const auto& section : surface.GetSections()) {
            for (const auto& design_var : section.GetDesignVariables()) {
                names.insert(design_var.name);
            }
        }
    }
    return std::vector<std::string>(names.begin(), names.end());
}

size_t Aircraft::GetNumDesignVariables() const {
    return GetDesignVariableNames().size();
}

int Aircraft::ComputeTotalVortexCount() const {
    int count = 0;
    for (const auto& surface : surfaces_) {
        int n_chord = surface.GetChordwisePanelCount();
        int n_span = surface.GetSpanwisePanelCount();
        count += n_chord * n_span;
        if (surface.IsYDuplicated()) {
            count *= 2;
        }
    }
    return count;
}

int Aircraft::ComputeTotalStripCount() const {
    int count = 0;
    for (const auto& surface : surfaces_) {
        count += surface.GetSpanwisePanelCount();
        if (surface.IsYDuplicated()) {
            count *= 2;
        }
    }
    return count;
}

Scalar Aircraft::ComputeTotalPlanformArea() const {
    Scalar area = 0.0;
    for (const auto& surface : surfaces_) {
        area += surface.ComputePlanformArea();
    }
    return area;
}

Scalar Aircraft::ComputeTotalSpan() const {
    Scalar max_span = 0.0;
    for (const auto& surface : surfaces_) {
        max_span = std::max(max_span, surface.ComputeSpan());
    }
    return max_span;
}

}  // namespace vlm
