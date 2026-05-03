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

#ifndef PEGASUS_GEOMETRY_BODY_H_
#define PEGASUS_GEOMETRY_BODY_H_

#include "pegasus/core/types.h"

#include <string>
#include <vector>

namespace pegasus {

/// @brief Cross-section definition for a body.
struct BodySection {
    Vec3 center;    ///< Center point of cross-section
    Scalar width;   ///< Half-width (Y extent)
    Scalar height;  ///< Half-height (Z extent)
};

/// @brief Fuselage body definition using source-doublet model.
///
/// Bodies in VLM are modeled using axial source and doublet distributions
/// along the centerline. This provides interference effects between the
/// fuselage and lifting surfaces.
///
/// A body is defined by a series of cross-sections along its length.
class Body {
public:
    Body() = default;

    /// @brief Creates a body with a name.
    explicit Body(std::string name);

    /// @name Basic Properties
    /// @{

    /// @brief Gets the body name.
    [[nodiscard]] const std::string& GetName() const { return name_; }

    /// @brief Sets the body name.
    void SetName(const std::string& name) { name_ = name; }

    /// @}

    /// @name Cross-Sections
    /// @{

    /// @brief Adds a cross-section.
    void AddSection(const BodySection& section) {
        sections_.push_back(section);
    }

    /// @brief Gets all cross-sections.
    [[nodiscard]] const std::vector<BodySection>& GetSections() const {
        return sections_;
    }

    /// @brief Gets the number of cross-sections.
    [[nodiscard]] size_t GetNumSections() const { return sections_.size(); }

    /// @}

    /// @name Paneling Parameters
    /// @{

    /// @brief Gets the number of source-doublet nodes.
    [[nodiscard]] int GetNodeCount() const { return node_count_; }

    /// @brief Sets the number of source-doublet nodes.
    void SetNodeCount(int count) { node_count_ = count; }

    /// @brief Gets the node spacing parameter.
    [[nodiscard]] Scalar GetNodeSpacing() const { return node_spacing_; }

    /// @brief Sets the node spacing parameter.
    void SetNodeSpacing(Scalar spacing) { node_spacing_ = spacing; }

    /// @}

    /// @name Coordinate Transformations
    /// @{

    /// @brief Gets the translation vector.
    [[nodiscard]] const Vec3& GetTranslation() const { return translation_; }

    /// @brief Sets the translation vector.
    void SetTranslation(const Vec3& translation) { translation_ = translation; }

    /// @brief Gets the scale factors.
    [[nodiscard]] const Vec3& GetScale() const { return scale_; }

    /// @brief Sets the scale factors.
    void SetScale(const Vec3& scale) { scale_ = scale; }

    /// @}

    /// @name Symmetry
    /// @{

    /// @brief Checks if body is Y-duplicated.
    [[nodiscard]] bool IsYDuplicated() const { return y_duplicate_; }

    /// @brief Sets Y-duplication flag.
    void SetYDuplicate(bool duplicate) { y_duplicate_ = duplicate; }

    /// @}

private:
    std::string name_;
    std::vector<BodySection> sections_;
    int node_count_{20};
    Scalar node_spacing_{1.0};
    Vec3 translation_{0, 0, 0};
    Vec3 scale_{1, 1, 1};
    bool y_duplicate_{false};
};

}  // namespace pegasus

#endif  // PEGASUS_GEOMETRY_BODY_H_
