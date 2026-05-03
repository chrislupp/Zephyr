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

#ifndef PEGASUS_GEOMETRY_AIRCRAFT_H_
#define PEGASUS_GEOMETRY_AIRCRAFT_H_

#include "pegasus/core/types.h"
#include "pegasus/geometry/body.h"
#include "pegasus/geometry/surface.h"

#include <memory>
#include <string>
#include <vector>

namespace pegasus {

/// @brief Complete aircraft configuration for VLM analysis.
///
/// An Aircraft holds the complete geometric definition including all
/// lifting surfaces, bodies, reference quantities, and symmetry settings.
/// This is the top-level geometry container loaded from AVL files.
class Aircraft {
public:
    Aircraft() = default;

    /// @brief Creates an aircraft with a title.
    explicit Aircraft(std::string title);

    /// @name Identification
    /// @{

    /// @brief Gets the aircraft title/name.
    [[nodiscard]] const std::string& GetTitle() const { return title_; }

    /// @brief Sets the aircraft title.
    void SetTitle(const std::string& title) { title_ = title; }

    /// @}

    /// @name Reference Quantities
    /// @{

    /// @brief Gets the reference area.
    [[nodiscard]] Scalar GetReferenceArea() const { return reference_area_; }

    /// @brief Sets the reference area.
    void SetReferenceArea(Scalar area) { reference_area_ = area; }

    /// @brief Gets the reference chord.
    [[nodiscard]] Scalar GetReferenceChord() const { return reference_chord_; }

    /// @brief Sets the reference chord.
    void SetReferenceChord(Scalar chord) { reference_chord_ = chord; }

    /// @brief Gets the reference span.
    [[nodiscard]] Scalar GetReferenceSpan() const { return reference_span_; }

    /// @brief Sets the reference span.
    void SetReferenceSpan(Scalar span) { reference_span_ = span; }

    /// @brief Gets the moment reference point.
    [[nodiscard]] const Vec3& GetMomentReferencePoint() const {
        return moment_reference_point_;
    }

    /// @brief Sets the moment reference point.
    void SetMomentReferencePoint(const Vec3& point) {
        moment_reference_point_ = point;
    }

    /// @}

    /// @name Flow Conditions
    /// @{

    /// @brief Gets the reference Mach number.
    [[nodiscard]] Scalar GetMachNumber() const { return mach_number_; }

    /// @brief Sets the reference Mach number.
    void SetMachNumber(Scalar mach) { mach_number_ = mach; }

    /// @brief Gets the reference profile drag coefficient.
    ///
    /// Constant drag coefficient added to induced drag.
    [[nodiscard]] Scalar GetProfileDragCoefficient() const {
        return cd_profile_;
    }

    /// @brief Sets the reference profile drag coefficient.
    void SetProfileDragCoefficient(Scalar cd) { cd_profile_ = cd; }

    /// @}

    /// @name Symmetry
    /// @{

    /// @brief Gets the Y-axis (XZ plane) symmetry type.
    [[nodiscard]] SymmetryType GetYSymmetry() const { return y_symmetry_; }

    /// @brief Sets the Y-axis symmetry type.
    void SetYSymmetry(SymmetryType sym) { y_symmetry_ = sym; }

    /// @brief Gets the Z-axis (XY plane) symmetry type.
    [[nodiscard]] SymmetryType GetZSymmetry() const { return z_symmetry_; }

    /// @brief Sets the Z-axis symmetry type.
    void SetZSymmetry(SymmetryType sym) { z_symmetry_ = sym; }

    /// @brief Gets the Z symmetry plane location.
    [[nodiscard]] Scalar GetZSymmetryPlane() const { return z_symmetry_plane_; }

    /// @brief Sets the Z symmetry plane location.
    void SetZSymmetryPlane(Scalar z) { z_symmetry_plane_ = z; }

    /// @}

    /// @name Surfaces
    /// @{

    /// @brief Adds a lifting surface.
    void AddSurface(const Surface& surface) { surfaces_.push_back(surface); }

    /// @brief Adds a lifting surface (move semantics).
    void AddSurface(Surface&& surface) {
        surfaces_.push_back(std::move(surface));
    }

    /// @brief Gets all lifting surfaces.
    [[nodiscard]] const std::vector<Surface>& GetSurfaces() const {
        return surfaces_;
    }

    /// @brief Gets mutable reference to surfaces.
    std::vector<Surface>& GetSurfacesMutable() { return surfaces_; }

    /// @brief Gets the number of surfaces.
    [[nodiscard]] size_t GetNumSurfaces() const { return surfaces_.size(); }

    /// @brief Gets a surface by index.
    [[nodiscard]] const Surface& GetSurface(size_t index) const {
        return surfaces_.at(index);
    }

    /// @brief Gets a surface by name.
    /// @return Pointer to surface, or nullptr if not found
    [[nodiscard]] const Surface* FindSurfaceByName(
        const std::string& name) const;

    /// @}

    /// @name Bodies
    /// @{

    /// @brief Adds a body.
    void AddBody(const Body& body) { bodies_.push_back(body); }

    /// @brief Gets all bodies.
    [[nodiscard]] const std::vector<Body>& GetBodies() const { return bodies_; }

    /// @brief Gets the number of bodies.
    [[nodiscard]] size_t GetNumBodies() const { return bodies_.size(); }

    /// @}

    /// @name Control Surfaces
    /// @{

    /// @brief Gets all control surface names across all surfaces.
    [[nodiscard]] std::vector<std::string> GetControlSurfaceNames() const;

    /// @brief Gets the number of unique control surfaces.
    [[nodiscard]] size_t GetNumControlSurfaces() const;

    /// @}

    /// @name Design Variables
    /// @{

    /// @brief Gets all design variable names across all surfaces.
    [[nodiscard]] std::vector<std::string> GetDesignVariableNames() const;

    /// @brief Gets the number of unique design variables.
    [[nodiscard]] size_t GetNumDesignVariables() const;

    /// @}

    /// @name Panel Counts
    /// @{

    /// @brief Computes the total number of vortex elements.
    [[nodiscard]] int ComputeTotalVortexCount() const;

    /// @brief Computes the total number of spanwise strips.
    [[nodiscard]] int ComputeTotalStripCount() const;

    /// @}

    /// @name Computed Properties
    /// @{

    /// @brief Computes total planform area from all surfaces.
    [[nodiscard]] Scalar ComputeTotalPlanformArea() const;

    /// @brief Computes total span from all surfaces.
    [[nodiscard]] Scalar ComputeTotalSpan() const;

    /// @}

private:
    std::string title_;

    // Reference quantities
    Scalar reference_area_{1.0};
    Scalar reference_chord_{1.0};
    Scalar reference_span_{1.0};
    Vec3 moment_reference_point_{0, 0, 0};

    // Flow conditions
    Scalar mach_number_{0.0};
    Scalar cd_profile_{0.0};

    // Symmetry
    SymmetryType y_symmetry_{SymmetryType::kNone};
    SymmetryType z_symmetry_{SymmetryType::kNone};
    Scalar z_symmetry_plane_{0.0};

    // Geometry
    std::vector<Surface> surfaces_;
    std::vector<Body> bodies_;
};

}  // namespace pegasus

#endif  // PEGASUS_GEOMETRY_AIRCRAFT_H_
