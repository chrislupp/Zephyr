// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_SOLVER_VORTEX_LATTICE_H_
#define VLM_SOLVER_VORTEX_LATTICE_H_

#include "vlm/core/types.h"
#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/mesh_input.h"
#include "vlm/solver/horseshoe_vortex.h"

#include <vector>

namespace vlm {

/// @brief Mesh indices for a surface in the vortex lattice.
struct SurfaceMeshInfo {
    /// Name of the surface
    std::string name;

    /// Index of first horseshoe vortex element in this surface
    Index first_vortex_index{0};

    /// Number of horseshoe vortex elements in this surface
    Index num_vortices{0};

    /// Index of first strip in this surface
    Index first_strip_index{0};

    /// Number of strips in this surface
    Index num_strips{0};

    /// Number of chordwise panels
    Index num_chordwise{0};

    /// Number of spanwise panels
    Index num_spanwise{0};

    /// Wetted area of surface (one side)
    Scalar wetted_area{0.0};

    /// Average chord of surface
    Scalar average_chord{0.0};

    /// Whether this is a y-duplicated (mirrored) surface
    bool is_y_duplicated{false};

    /// Index of original surface (if y-duplicated)
    Index original_surface_index{0};
};

/// @brief Contains the complete vortex lattice mesh for an aircraft.
///
/// The VortexLattice class discretizes an Aircraft geometry into a mesh of
/// horseshoe vortex elements organized into spanwise strips. It provides
/// the geometric data needed for the VLM solver.
///
/// The mesh generation follows the AVL methodology:
/// - Sections define spanwise stations with chord, twist, and airfoil
/// - Cosine spacing distributes panels with finer resolution at edges
/// - Y-duplicate surfaces are reflected about a symmetry plane
/// - Control surfaces are mapped to element influence fractions
class VortexLattice {
public:
    VortexLattice() = default;

    /// @brief Generates the vortex lattice mesh from an aircraft configuration.
    ///
    /// @param aircraft The aircraft geometry to discretize
    /// @throws std::runtime_error if mesh generation fails
    void GenerateFromAircraft(const Aircraft& aircraft);

    /// @brief Generates the vortex lattice mesh from a direct mesh input.
    ///
    /// Each SurfaceMesh in the input defines a structured grid of node
    /// coordinates. Panels are formed from adjacent nodes, with bound
    /// vortices at 1/4-chord and control points at 3/4-chord.
    ///
    /// @param input The mesh input containing surface node grids
    void GenerateFromMesh(const MeshInput& input);

    /// @brief Gets the stored mesh input (valid only after GenerateFromMesh).
    const MeshInput* GetMeshInput() const {
        return has_mesh_input_ ? &mesh_input_ : nullptr;
    }

    /// @brief Gets the total number of horseshoe vortex elements.
    Index GetNumVortices() const {
        return static_cast<Index>(vortices_.size());
    }

    /// @brief Gets the total number of spanwise strips.
    Index GetNumStrips() const { return static_cast<Index>(strips_.size()); }

    /// @brief Gets the number of surfaces in the mesh.
    Index GetNumSurfaces() const {
        return static_cast<Index>(surface_info_.size());
    }

    /// @brief Gets a horseshoe vortex element by index.
    const HorseshoeVortex& GetVortex(Index index) const {
        return vortices_[index];
    }
    HorseshoeVortex& GetVortex(Index index) { return vortices_[index]; }

    /// @brief Gets a strip by index.
    const Strip& GetStrip(Index index) const { return strips_[index]; }
    Strip& GetStrip(Index index) { return strips_[index]; }

    /// @brief Gets mesh information for a surface.
    const SurfaceMeshInfo& GetSurfaceInfo(Index index) const {
        return surface_info_[index];
    }

    /// @brief Gets all horseshoe vortex elements.
    const std::vector<HorseshoeVortex>& GetVortices() const {
        return vortices_;
    }
    std::vector<HorseshoeVortex>& GetVortices() { return vortices_; }

    /// @brief Gets all strips.
    const std::vector<Strip>& GetStrips() const { return strips_; }
    std::vector<Strip>& GetStrips() { return strips_; }

    /// @brief Gets control surface influences for all elements.
    const std::vector<std::vector<ControlSurfaceInfluence>>&
    GetControlSurfaceInfluences() const {
        return control_influences_;
    }

    /// @brief Gets control surface names defined in the mesh.
    const std::vector<std::string>& GetControlSurfaceNames() const {
        return control_surface_names_;
    }

    /// @brief Clears all mesh data.
    void Clear();

    /// @brief Recalculates normal vectors (after geometry changes).
    void RecalculateNormals();

    /// @brief Gets the vortex core radius for control point calculations.
    Scalar GetVortexCoreRadiusChordFraction() const {
        return vortex_core_radius_chord_fraction_;
    }

    /// @brief Sets the vortex core radius as fraction of chord.
    void SetVortexCoreRadiusChordFraction(Scalar fraction) {
        vortex_core_radius_chord_fraction_ = fraction;
    }

    /// @brief Gets the vortex core radius for wake calculations.
    Scalar GetVortexCoreRadiusWidthFraction() const {
        return vortex_core_radius_width_fraction_;
    }

    /// @brief Sets the vortex core radius as fraction of strip width.
    void SetVortexCoreRadiusWidthFraction(Scalar fraction) {
        vortex_core_radius_width_fraction_ = fraction;
    }

private:
    /// Generates mesh for a single surface from Aircraft sections.
    void GenerateSurfaceMesh(const Surface& surface, Index surface_index,
                             const Aircraft& aircraft);

    /// Generates mesh for a single surface from a structured mesh grid.
    void GenerateSurfaceFromMesh(const SurfaceMesh& surface_mesh,
                                 Index surface_index);

    /// Generates a y-duplicated (mirrored) surface mesh.
    void GenerateDuplicateSurfaceMesh(Index original_surface_index,
                                      Scalar y_duplicate_plane);

    /// Calculates normal vectors for all elements.
    void CalculateNormals();

    /// Calculates chordwise spacing fractions for vortex and control points.
    void CalculateChordwiseSpacing(Index num_chordwise,
                                   SpacingType spacing_type,
                                   Scalar cl_alpha_factor,
                                   std::vector<Scalar>& x_vortex,
                                   std::vector<Scalar>& x_control) const;

    /// Horseshoe vortex elements
    std::vector<HorseshoeVortex> vortices_;

    /// Spanwise strips
    std::vector<Strip> strips_;

    /// Mesh information for each surface
    std::vector<SurfaceMeshInfo> surface_info_;

    /// Control surface influences for each vortex element
    std::vector<std::vector<ControlSurfaceInfluence>> control_influences_;

    /// Names of control surfaces in the mesh
    std::vector<std::string> control_surface_names_;

    /// Vortex core radius as fraction of local chord (for control point calcs)
    Scalar vortex_core_radius_chord_fraction_{0.75};

    /// Vortex core radius as fraction of strip width (for wake calcs)
    Scalar vortex_core_radius_width_fraction_{0.25};

    /// Spanwise axis fraction for perpendicular velocity definition
    Scalar spanwise_axis_fraction_{0.25};

    /// Stored mesh input for gradient computation
    MeshInput mesh_input_;
    bool has_mesh_input_{false};
};

}  // namespace vlm

#endif  // VLM_SOLVER_VORTEX_LATTICE_H_
