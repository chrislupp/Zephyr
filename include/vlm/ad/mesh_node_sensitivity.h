// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#pragma once

/// @file mesh_node_sensitivity.h
/// @brief Chain rule from VLM mesh positions to structured mesh node coordinates.
///
/// When using direct mesh input (MeshInput), the mesh node coordinates are
/// the design variables. This module provides the reverse-mode chain rule
/// from MeshGradients (d/d rv1, rv2, rc, normal) back to mesh node coordinates.
///
/// Forward relationships (for panel (i,j)):
///   rv1 = 0.75*N(i,j) + 0.25*N(i+1,j)          [left 1/4 chord]
///   rv2 = 0.75*N(i,j+1) + 0.25*N(i+1,j+1)      [right 1/4 chord]
///   rc  = 0.5*(0.25*N(i,j) + 0.75*N(i+1,j))
///       + 0.5*(0.25*N(i,j+1) + 0.75*N(i+1,j+1)) [midspan 3/4 chord]
///   normal = normalize((N(i+1,j+1)-N(i,j)) x (N(i+1,j)-N(i,j+1)))

#include <cmath>
#include <vector>

#include "vlm/ad/mesh_sensitivity.h"
#include "vlm/core/types.h"
#include "vlm/geometry/mesh_input.h"

namespace vlm {
namespace ad {

/// @brief Gradients of aerodynamic outputs w.r.t. mesh node coordinates.
///
/// For a surface mesh with nx chordwise and ny spanwise panels,
/// the gradient array has size (nx+1)*(ny+1)*3 corresponding to
/// the x,y,z coordinates of each node.
struct MeshNodeGradients {
  Index nx = 0;
  Index ny = 0;
  VecX d_nodes;  ///< [(nx+1)*(ny+1)*3] flat array

  MeshNodeGradients() = default;
  MeshNodeGradients(Index nx_, Index ny_)
      : nx(nx_), ny(ny_), d_nodes(VecX::Zero((nx_ + 1) * (ny_ + 1) * 3)) {}

  /// Access d/d(Node(i,j).k) where k=0,1,2 for x,y,z
  double& operator()(Index i, Index j, Index k) {
    return d_nodes[(i * (ny + 1) + j) * 3 + k];
  }

  double operator()(Index i, Index j, Index k) const {
    return d_nodes[(i * (ny + 1) + j) * 3 + k];
  }

  /// Add contribution to node (i,j) from a 3-vector
  void AddToNode(Index i, Index j, const double* d, double scale = 1.0) {
    Index base = (i * (ny + 1) + j) * 3;
    d_nodes[base + 0] += scale * d[0];
    d_nodes[base + 1] += scale * d[1];
    d_nodes[base + 2] += scale * d[2];
  }

  void Zero() { d_nodes.setZero(); }

  Index NumNodes() const { return (nx + 1) * (ny + 1); }
};

/// @brief All mesh node gradients for all aero outputs on one surface.
struct SurfaceMeshNodeGradients {
  std::string surface_name;
  Index nx = 0;
  Index ny = 0;

  MeshNodeGradients dCL;
  MeshNodeGradients dCD;
  MeshNodeGradients dCY;
  MeshNodeGradients dCl;
  MeshNodeGradients dCm;
  MeshNodeGradients dCn;

  void Init(const std::string& name, Index nx_, Index ny_) {
    surface_name = name;
    nx = nx_;
    ny = ny_;
    dCL = MeshNodeGradients(nx_, ny_);
    dCD = MeshNodeGradients(nx_, ny_);
    dCY = MeshNodeGradients(nx_, ny_);
    dCl = MeshNodeGradients(nx_, ny_);
    dCm = MeshNodeGradients(nx_, ny_);
    dCn = MeshNodeGradients(nx_, ny_);
  }
};

// ============================================================================
// Normal vector derivative kernel
// ============================================================================

/// @brief Propagate d_normal gradient back through normalize(d1 x d2) to
/// the four corner nodes of a panel.
///
/// Forward:
///   d1 = N(i+1,j+1) - N(i,j)
///   d2 = N(i+1,j) - N(i,j+1)
///   raw = d1 x d2
///   normal = raw / |raw|
///
/// Adjoint: given d_normal[3], compute contributions to the 4 corner nodes.
///
/// @param d_normal Input: d(objective)/d(normal) [3]
/// @param n00 Node(i,j) coordinates [3]
/// @param n10 Node(i+1,j) coordinates [3]
/// @param n01 Node(i,j+1) coordinates [3]
/// @param n11 Node(i+1,j+1) coordinates [3]
/// @param d_n00 Output: accumulated gradient for Node(i,j) [3]
/// @param d_n10 Output: accumulated gradient for Node(i+1,j) [3]
/// @param d_n01 Output: accumulated gradient for Node(i,j+1) [3]
/// @param d_n11 Output: accumulated gradient for Node(i+1,j+1) [3]
inline void PropagateNormalGradientToNodes(
    const double* d_normal,
    const double* n00, const double* n10,
    const double* n01, const double* n11,
    double* d_n00, double* d_n10,
    double* d_n01, double* d_n11) {
  // Forward: d1 = n11 - n00, d2 = n10 - n01
  double d1[3] = {n11[0] - n00[0], n11[1] - n00[1], n11[2] - n00[2]};
  double d2[3] = {n10[0] - n01[0], n10[1] - n01[1], n10[2] - n01[2]};

  // raw = d1 x d2
  double raw[3];
  raw[0] = d1[1] * d2[2] - d1[2] * d2[1];
  raw[1] = d1[2] * d2[0] - d1[0] * d2[2];
  raw[2] = d1[0] * d2[1] - d1[1] * d2[0];

  double mag = std::sqrt(raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]);
  if (mag < 1e-12) return;

  double inv_mag = 1.0 / mag;
  double normal[3] = {raw[0] * inv_mag, raw[1] * inv_mag, raw[2] * inv_mag};

  // Adjoint of normalization: d_raw = (d_normal - normal*(normal.d_normal)) / mag
  double n_dot_dn = normal[0] * d_normal[0] + normal[1] * d_normal[1] +
                    normal[2] * d_normal[2];
  double d_raw[3];
  d_raw[0] = (d_normal[0] - normal[0] * n_dot_dn) * inv_mag;
  d_raw[1] = (d_normal[1] - normal[1] * n_dot_dn) * inv_mag;
  d_raw[2] = (d_normal[2] - normal[2] * n_dot_dn) * inv_mag;

  // Adjoint of cross product: raw = d1 x d2
  // d_d1 += d_raw x d2  (actually d2 x d_raw for adjoint, but cross product adjoint is:
  //   d_d1 += d2 x d_raw  ... wait, let me be careful.
  //   raw = d1 x d2
  //   d(raw)/d(d1) applied to d_raw gives: d_d1 += d_raw x d2 ... no.
  //   Actually: d(d1 x d2)/d(d1) . delta_d1 = delta_d1 x d2
  //   So adjoint: d_d1 += -d2 x d_raw = d_raw x d2 ... hmm.
  //
  // Let's use components. raw[0] = d1[1]*d2[2] - d1[2]*d2[1]
  // d_d1[1] += d2[2] * d_raw[0], d_d1[2] += -d2[1] * d_raw[0]
  // etc. This is equivalent to: d_d1 += cross(d2, d_raw) ... no.
  // d_d1 from raw = d1 x d2:  d_d1 = -d2 x d_raw  (i.e., d2 crossed with d_raw, negated)
  // Actually: d_d1[k] = sum_m d(raw_m)/d(d1_k) * d_raw[m]
  // raw[0] = d1[1]*d2[2] - d1[2]*d2[1]  => d/d(d1[1]) = d2[2], d/d(d1[2]) = -d2[1]
  // raw[1] = d1[2]*d2[0] - d1[0]*d2[2]  => d/d(d1[2]) = d2[0], d/d(d1[0]) = -d2[2]
  // raw[2] = d1[0]*d2[1] - d1[1]*d2[0]  => d/d(d1[0]) = d2[1], d/d(d1[1]) = -d2[0]
  //
  // d_d1[0] = -d2[2]*d_raw[1] + d2[1]*d_raw[2]  = (d2 x d_raw)[0]... wait
  // d2 x d_raw = (d2[1]*d_raw[2]-d2[2]*d_raw[1], d2[2]*d_raw[0]-d2[0]*d_raw[2], d2[0]*d_raw[1]-d2[1]*d_raw[0])
  // So d_d1[0] = d2[1]*d_raw[2] - d2[2]*d_raw[1] ... but we had d_d1[0] = -d2[2]*d_raw[1] + d2[1]*d_raw[2]
  // Yes! So d_d1 = d2 x d_raw. Good, and symmetrically d_d2 = d_raw x d1.

  double d_d1[3], d_d2[3];
  // d_d1 = d2 x d_raw
  d_d1[0] = d2[1] * d_raw[2] - d2[2] * d_raw[1];
  d_d1[1] = d2[2] * d_raw[0] - d2[0] * d_raw[2];
  d_d1[2] = d2[0] * d_raw[1] - d2[1] * d_raw[0];

  // d_d2 = d_raw x d1
  d_d2[0] = d_raw[1] * d1[2] - d_raw[2] * d1[1];
  d_d2[1] = d_raw[2] * d1[0] - d_raw[0] * d1[2];
  d_d2[2] = d_raw[0] * d1[1] - d_raw[1] * d1[0];

  // Adjoint of d1 = n11 - n00: d_n11 += d_d1, d_n00 -= d_d1
  // Adjoint of d2 = n10 - n01: d_n10 += d_d2, d_n01 -= d_d2
  for (int k = 0; k < 3; ++k) {
    d_n11[k] += d_d1[k];
    d_n00[k] -= d_d1[k];
    d_n10[k] += d_d2[k];
    d_n01[k] -= d_d2[k];
  }
}

// ============================================================================
// Linear position chain rule (scatter-add)
// ============================================================================

/// @brief Propagate MeshGradients to mesh node gradients for one surface.
///
/// For each panel (i,j), the vortex positions are linear combinations of
/// the 4 corner nodes. The normal is nonlinear (cross product + normalize).
///
/// @param mesh_grads MeshGradients (d_rv1, d_rv2, d_rc, d_normal)
/// @param surface_mesh Original surface mesh (for node coordinates)
/// @param first_vortex First vortex index for this surface
/// @param nx Chordwise panels
/// @param ny Spanwise panels
/// @param out Output node gradients (accumulated)
inline void PropagateMeshGradientsToNodes(
    const MeshGradients& mesh_grads,
    const SurfaceMesh& surface_mesh,
    Index first_vortex,
    Index nx, Index ny,
    MeshNodeGradients& out) {
  for (Index j = 0; j < ny; ++j) {
    for (Index i = 0; i < nx; ++i) {
      Index v = first_vortex + j * nx + i;

      const double* d_rv1 = &mesh_grads.d_rv1[3 * v];
      const double* d_rv2 = &mesh_grads.d_rv2[3 * v];
      const double* d_rc = &mesh_grads.d_rc[3 * v];
      const double* d_normal = &mesh_grads.d_normal[3 * v];

      // rv1 = 0.75*N(i,j) + 0.25*N(i+1,j)
      out.AddToNode(i, j, d_rv1, 0.75);
      out.AddToNode(i + 1, j, d_rv1, 0.25);

      // rv2 = 0.75*N(i,j+1) + 0.25*N(i+1,j+1)
      out.AddToNode(i, j + 1, d_rv2, 0.75);
      out.AddToNode(i + 1, j + 1, d_rv2, 0.25);

      // rc = 0.5*(0.25*N(i,j) + 0.75*N(i+1,j)) + 0.5*(0.25*N(i,j+1) + 0.75*N(i+1,j+1))
      out.AddToNode(i, j, d_rc, 0.125);       // 0.5 * 0.25
      out.AddToNode(i + 1, j, d_rc, 0.375);   // 0.5 * 0.75
      out.AddToNode(i, j + 1, d_rc, 0.125);
      out.AddToNode(i + 1, j + 1, d_rc, 0.375);

      // Normal: nonlinear — propagate through cross product + normalize
      double n00[3], n10[3], n01[3], n11[3];
      Vec3 p00 = surface_mesh.GetNode(i, j);
      Vec3 p10 = surface_mesh.GetNode(i + 1, j);
      Vec3 p01 = surface_mesh.GetNode(i, j + 1);
      Vec3 p11 = surface_mesh.GetNode(i + 1, j + 1);
      for (int k = 0; k < 3; ++k) {
        n00[k] = p00[k];
        n10[k] = p10[k];
        n01[k] = p01[k];
        n11[k] = p11[k];
      }

      double d_n00[3] = {0, 0, 0};
      double d_n10[3] = {0, 0, 0};
      double d_n01[3] = {0, 0, 0};
      double d_n11[3] = {0, 0, 0};

      PropagateNormalGradientToNodes(d_normal, n00, n10, n01, n11,
                                     d_n00, d_n10, d_n01, d_n11);

      out.AddToNode(i, j, d_n00);
      out.AddToNode(i + 1, j, d_n10);
      out.AddToNode(i, j + 1, d_n01);
      out.AddToNode(i + 1, j + 1, d_n11);
    }
  }
}

}  // namespace ad
}  // namespace vlm
