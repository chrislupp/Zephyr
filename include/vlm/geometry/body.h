// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_GEOMETRY_BODY_H_
#define VLM_GEOMETRY_BODY_H_

#include "vlm/core/types.h"

#include <string>
#include <vector>

namespace vlm {

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
    explicit Body(const std::string& name);

    /// @name Basic Properties
    /// @{

    /// @brief Gets the body name.
    const std::string& GetName() const { return name_; }

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
    const std::vector<BodySection>& GetSections() const { return sections_; }

    /// @brief Gets the number of cross-sections.
    size_t GetNumSections() const { return sections_.size(); }

    /// @}

    /// @name Paneling Parameters
    /// @{

    /// @brief Gets the number of source-doublet nodes.
    int GetNodeCount() const { return node_count_; }

    /// @brief Sets the number of source-doublet nodes.
    void SetNodeCount(int count) { node_count_ = count; }

    /// @brief Gets the node spacing parameter.
    Scalar GetNodeSpacing() const { return node_spacing_; }

    /// @brief Sets the node spacing parameter.
    void SetNodeSpacing(Scalar spacing) { node_spacing_ = spacing; }

    /// @}

    /// @name Coordinate Transformations
    /// @{

    /// @brief Gets the translation vector.
    const Vec3& GetTranslation() const { return translation_; }

    /// @brief Sets the translation vector.
    void SetTranslation(const Vec3& translation) { translation_ = translation; }

    /// @brief Gets the scale factors.
    const Vec3& GetScale() const { return scale_; }

    /// @brief Sets the scale factors.
    void SetScale(const Vec3& scale) { scale_ = scale; }

    /// @}

    /// @name Symmetry
    /// @{

    /// @brief Checks if body is Y-duplicated.
    bool IsYDuplicated() const { return y_duplicate_; }

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

}  // namespace vlm

#endif  // VLM_GEOMETRY_BODY_H_
