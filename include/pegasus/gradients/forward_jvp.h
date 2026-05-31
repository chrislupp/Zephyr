/*
    Project Pegasus
    Version 0.2.0

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

#pragma once

/// @file forward_jvp.h
/// @brief Forward-mode Jacobian-vector products (JVPs) for the VLM system.
///
/// Provides forward-mode derivatives through the VLM solve chain for use in
/// coupled aerostructural condensed Jacobian computation. The chain is:
///
///   mesh_nodes -> vortex_geometry -> (AIC, b) -> Gamma -> forces
///
/// Given a perturbation in mesh nodes (delta_nodes), these functions compute:
///   1. Vortex geometry perturbation (delta_rv1, delta_rv2, delta_rc,
///   delta_normal)
///   2. VLM residual perturbation: delta(AIC * Gamma - b)
///   3. Total per-vortex force perturbation (including implicit VLM solve)

#include "pegasus/core/types.h"
#include "pegasus/geometry/mesh_input.h"
#include "pegasus/gradients/mesh_node_sensitivity.h"
#include "pegasus/solver/vortex_lattice_solver.h"

namespace pegasus::gradients {

/// @brief Perturbation of vortex geometry positions.
///
/// Stores forward-mode perturbations of the four geometric quantities that
/// define VLM panel geometry. These are the forward-mode counterpart of
/// MeshGradients (which stores reverse-mode adjoint seeds).
struct VortexGeometryPerturbation {
    VecX d_rv1;     ///< Left bound vortex endpoints [3*N]
    VecX d_rv2;     ///< Right bound vortex endpoints [3*N]
    VecX d_rc;      ///< Control points [3*N]
    VecX d_normal;  ///< Panel normals [3*N]

    void Resize(Index n) {
        d_rv1 = VecX::Zero(3 * n);
        d_rv2 = VecX::Zero(3 * n);
        d_rc = VecX::Zero(3 * n);
        d_normal = VecX::Zero(3 * n);
    }
};

/// @brief Forward propagate mesh node perturbations to vortex geometry.
///
/// This is the forward mode of PropagateMeshGradientsToNodes. For each panel
/// (i,j) with vortex index v = first_vortex + j*nx + i:
///   delta_rv1[v] = 0.75*dN(i,j) + 0.25*dN(i+1,j)
///   delta_rv2[v] = 0.75*dN(i,j+1) + 0.25*dN(i+1,j+1)
///   delta_rc[v]  = 0.125*dN(i,j) + 0.375*dN(i+1,j) + 0.125*dN(i,j+1) +
///                  0.375*dN(i+1,j+1)
///   delta_normal[v] = forward through normalize(cross(d1, d2))
///
/// @param surface_mesh Surface mesh (for node coordinates in normal derivative)
/// @param first_vortex First vortex index for this surface
/// @param nx Chordwise panels
/// @param ny Spanwise panels
/// @param delta_nodes Perturbation of mesh node coordinates
/// @return Vortex geometry perturbation
VortexGeometryPerturbation PropagateNodePerturbationsToVortex(
    const SurfaceMesh& surface_mesh, Index first_vortex, Index nx, Index ny,
    const MeshNodeGradients& delta_nodes);

/// @brief Compute forward-mode JVP of VLM residual w.r.t. vortex geometry.
///
/// Computes delta(AIC * Gamma - b) given perturbations in vortex geometry.
///
/// For each row i of the VLM system:
///   delta_R[i] = sum_j Gamma_j * delta(A[i,j]) - V_inf . delta_normal[i]
///
/// where delta(A[i,j]) accounts for perturbations in rc_i, rv1_j, rv2_j,
/// and normal_i using analytical Biot-Savart derivatives.
///
/// @param solver VLM solver (provides mesh, AIC parameters, and LU
/// factorization)
/// @param solution Forward solution (circulation values)
/// @param condition Flight condition (for freestream direction)
/// @param delta_geom Perturbation of vortex geometry
/// @return Perturbation of VLM residual (N x 1)
VecX ComputeResidualGeometryJVP(const VortexLatticeSolver& solver,
                                const SolutionState& solution,
                                const FlightCondition& condition,
                                const VortexGeometryPerturbation& delta_geom);

/// @brief Compute total per-vortex force JVP including VLM chain.
///
/// Computes the total perturbation of per-vortex body-frame forces given
/// vortex geometry perturbation and circulation perturbation:
///
///   delta_F_i = 2*Gamma_i*(delta_V_eff_i x dl_i + V_eff_i x delta_dl_i)
///             + 2*delta_Gamma_i*(V_eff_i x dl_i)
///
/// where:
///   delta_dl_i = delta_rv2_i - delta_rv1_i  (from delta_geom)
///   delta_V_eff_i = delta_VIM_i * Gamma + VIM_i * delta_Gamma
///   delta_VIM comes from vortex geometry perturbation (Biot-Savart
///   derivatives)
///
/// @param solver VLM solver (mesh, VIM)
/// @param solution Forward solution (circulation, induced velocity)
/// @param condition Flight condition
/// @param delta_geom Vortex geometry perturbation
/// @param delta_gamma Circulation perturbation (N x 1, from forward AIC solve)
/// @return Per-vortex force perturbation (3N x 1)
VecX ComputePerVortexForceTotalJVP(const VortexLatticeSolver& solver,
                                   const SolutionState& solution,
                                   const FlightCondition& condition,
                                   const VortexGeometryPerturbation& delta_geom,
                                   const VecX& delta_gamma);

}  // namespace pegasus::gradients
