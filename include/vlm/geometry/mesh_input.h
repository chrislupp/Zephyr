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

#ifndef VLM_GEOMETRY_MESH_INPUT_H_
#define VLM_GEOMETRY_MESH_INPUT_H_

#include "vlm/core/types.h"

#include <string>
#include <vector>

namespace vlm {

/// @brief Reference quantities for aerodynamic coefficient normalization.
///
/// These are needed by the solver and force calculator when using
/// direct mesh input (where no Aircraft object is available).
struct ReferenceQuantities {
    Scalar reference_area = 1.0;
    Scalar reference_chord = 1.0;
    Scalar reference_span = 1.0;
    Vec3 moment_reference_point = {0, 0, 0};
};

/// @brief A structured surface mesh for direct mesh input.
///
/// The mesh is a structured grid of (nx+1) x (ny+1) nodes, where nx is
/// the number of chordwise panels and ny is the number of spanwise panels.
/// Nodes are stored row-major: Node(i,j) is at row i*(ny+1)+j in the
/// nodes matrix. i=0 is the leading edge, i=nx is the trailing edge.
///
/// This follows the OpenAeroStruct convention where the mesh defines
/// the wing surface geometry directly, rather than using section parameters
/// (chord, twist, etc.).
struct SurfaceMesh {
    /// Name of the surface
    std::string name;

    /// Number of chordwise panels (mesh has nx+1 chordwise node rows)
    Index nx = 0;

    /// Number of spanwise panels (mesh has ny+1 spanwise node columns)
    Index ny = 0;

    /// Node coordinates: (nx+1)*(ny+1) rows x 3 columns
    /// Node(i,j) is at row i*(ny+1)+j. i=0 is LE, i=nx is TE.
    MatX nodes;

    /// Component index for vortex core radius logic
    Index component_index = 0;

    /// Whether elements respond to freestream alpha/beta
    bool affected_by_freestream = true;

    /// Whether this surface sheds a wake
    bool has_wake = true;

    /// Whether to create a y-duplicated (mirrored) copy
    bool y_duplicated = false;

    /// Y-plane for duplication (typically 0)
    Scalar y_duplicate_plane = 0.0;

    /// Access node (i, j) as a 3D point
    Vec3 GetNode(Index i, Index j) const {
        Index row = i * (ny + 1) + j;
        return nodes.row(row).transpose();
    }

    /// Set node (i, j)
    void SetNode(Index i, Index j, const Vec3& point) {
        Index row = i * (ny + 1) + j;
        nodes.row(row) = point.transpose();
    }

    /// Total number of nodes
    Index NumNodes() const { return (nx + 1) * (ny + 1); }

    /// Total number of panels
    Index NumPanels() const { return nx * ny; }
};

/// @brief Complete mesh-based input for VLM analysis.
///
/// This is the alternative to Aircraft for cases where the user wants
/// to provide mesh node coordinates directly rather than section-based
/// geometry parameters (chord, twist, LE position).
struct MeshInput {
    /// Surface meshes
    std::vector<SurfaceMesh> surfaces;

    /// Reference quantities for coefficient normalization
    ReferenceQuantities reference;
};

}  // namespace vlm

#endif  // VLM_GEOMETRY_MESH_INPUT_H_
