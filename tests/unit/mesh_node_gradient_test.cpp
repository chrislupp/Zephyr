// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include <gtest/gtest.h>

#include <cmath>

#include "vlm/ad/gradient_calculator.h"
#include "vlm/ad/mesh_node_sensitivity.h"
#include "vlm/aero/force_calculator.h"
#include "vlm/geometry/mesh_input.h"
#include "vlm/solver/vortex_lattice_solver.h"

namespace vlm {
namespace ad {
namespace {

/// Build a flat rectangular wing mesh.
SurfaceMesh BuildRectangularMesh(Index nx, Index ny,
                                  Scalar chord, Scalar span) {
  SurfaceMesh sm;
  sm.name = "Wing";
  sm.nx = nx;
  sm.ny = ny;
  sm.nodes.resize((nx + 1) * (ny + 1), 3);

  for (Index i = 0; i <= nx; ++i) {
    for (Index j = 0; j <= ny; ++j) {
      Scalar x = chord * static_cast<Scalar>(i) / static_cast<Scalar>(nx);
      Scalar y = span * static_cast<Scalar>(j) / static_cast<Scalar>(ny);
      Scalar z = 0.0;
      sm.SetNode(i, j, Vec3(x, y, z));
    }
  }
  return sm;
}

MeshInput BuildTestInput(Index nx, Index ny, Scalar chord, Scalar span) {
  SurfaceMesh sm = BuildRectangularMesh(nx, ny, chord, span);
  MeshInput input;
  input.surfaces.push_back(sm);
  input.reference.reference_area = chord * span;
  input.reference.reference_chord = chord;
  input.reference.reference_span = span;
  input.reference.moment_reference_point = Vec3(0.25 * chord, 0.5 * span, 0);
  return input;
}

/// Solve the VLM problem and return aero results
AerodynamicResults SolveAndGetResults(const MeshInput& input,
                                       const FlightCondition& condition) {
  VortexLatticeSolver solver;
  solver.SetMesh(input);
  SolutionState state = solver.Solve(condition);
  ForceCalculator calc;
  return calc.ComputeForces(input.reference, solver.GetMesh(), state, condition);
}

// Test that ComputeMeshNodeGradients returns correct structure
TEST(MeshNodeGradientTest, OutputStructure) {
  const Index nx = 2;
  const Index ny = 4;
  MeshInput input = BuildTestInput(nx, ny, 1.0, 4.0);

  GradientCalculator calc;
  calc.SetMesh(input);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  auto grads = calc.ComputeMeshNodeGradients(condition);

  ASSERT_EQ(grads.size(), 1u);  // One surface
  EXPECT_EQ(grads[0].nx, nx);
  EXPECT_EQ(grads[0].ny, ny);
  EXPECT_EQ(grads[0].dCL.d_nodes.size(), (nx + 1) * (ny + 1) * 3);
  EXPECT_EQ(grads[0].dCD.d_nodes.size(), (nx + 1) * (ny + 1) * 3);
}

// Finite difference validation of CL gradient w.r.t. mesh nodes
TEST(MeshNodeGradientTest, CLGradientFiniteDifference) {
  const Index nx = 2;
  const Index ny = 3;
  const Scalar chord = 1.0;
  const Scalar span = 3.0;
  MeshInput input = BuildTestInput(nx, ny, chord, span);

  FlightCondition condition;
  condition.alpha = 5.0 * kDegreesToRadians;

  // Compute analytical gradients
  GradientCalculator calc;
  calc.SetMesh(input);
  auto grads = calc.ComputeMeshNodeGradients(condition);
  ASSERT_EQ(grads.size(), 1u);

  const MeshNodeGradients& dCL = grads[0].dCL;
  const MeshNodeGradients& dCD = grads[0].dCD;

  // Finite difference validation
  const Scalar h = 1e-6;
  Index num_nodes = (nx + 1) * (ny + 1);
  int num_validated = 0;
  int num_failed = 0;

  for (Index node = 0; node < num_nodes; ++node) {
    for (int dir = 0; dir < 3; ++dir) {
      // Perturb +h
      MeshInput input_plus = input;
      input_plus.surfaces[0].nodes(node, dir) += h;
      AerodynamicResults res_plus = SolveAndGetResults(input_plus, condition);

      // Perturb -h
      MeshInput input_minus = input;
      input_minus.surfaces[0].nodes(node, dir) -= h;
      AerodynamicResults res_minus = SolveAndGetResults(input_minus, condition);

      Scalar fd_dCL = (res_plus.CL - res_minus.CL) / (2.0 * h);
      Scalar fd_dCD = (res_plus.CD - res_minus.CD) / (2.0 * h);

      Scalar anal_dCL = dCL.d_nodes[node * 3 + dir];
      Scalar anal_dCD = dCD.d_nodes[node * 3 + dir];

      // Use relative tolerance for larger values, absolute for small ones
      Scalar tol_CL = std::max(1e-4, 1e-3 * std::abs(fd_dCL));
      Scalar tol_CD = std::max(1e-4, 1e-3 * std::abs(fd_dCD));

      if (std::abs(anal_dCL - fd_dCL) > tol_CL) {
        num_failed++;
        if (num_failed <= 5) {
          EXPECT_NEAR(anal_dCL, fd_dCL, tol_CL)
              << "CL gradient mismatch at node=" << node << " dir=" << dir
              << " analytical=" << anal_dCL << " fd=" << fd_dCL;
        }
      }
      if (std::abs(anal_dCD - fd_dCD) > tol_CD) {
        num_failed++;
        if (num_failed <= 5) {
          EXPECT_NEAR(anal_dCD, fd_dCD, tol_CD)
              << "CD gradient mismatch at node=" << node << " dir=" << dir
              << " analytical=" << anal_dCD << " fd=" << fd_dCD;
        }
      }
      num_validated += 2;
    }
  }

  // At least most gradients should pass
  double pass_rate = static_cast<double>(num_validated - num_failed) /
                     static_cast<double>(num_validated);
  EXPECT_GT(pass_rate, 0.8) << "Only " << (pass_rate * 100)
                             << "% of gradients passed FD validation ("
                             << num_failed << "/" << num_validated << " failed)";
}

// Test normal gradient propagation kernel
TEST(MeshNodeGradientTest, NormalGradientKernel) {
  // Create a simple non-degenerate panel
  double n00[3] = {0, 0, 0};
  double n10[3] = {1, 0, 0};
  double n01[3] = {0, 1, 0};
  double n11[3] = {1, 1, 0.1};  // Slight camber

  // Seed: d_normal = (0, 0, 1)
  double d_normal[3] = {0, 0, 1};

  double d_n00[3] = {0, 0, 0};
  double d_n10[3] = {0, 0, 0};
  double d_n01[3] = {0, 0, 0};
  double d_n11[3] = {0, 0, 0};

  PropagateNormalGradientToNodes(d_normal, n00, n10, n01, n11,
                                 d_n00, d_n10, d_n01, d_n11);

  // Verify with finite differences
  const double h = 1e-7;
  auto compute_normal = [](const double* p00, const double* p10,
                           const double* p01, const double* p11) -> Vec3 {
    Vec3 d1(p11[0]-p00[0], p11[1]-p00[1], p11[2]-p00[2]);
    Vec3 d2(p10[0]-p01[0], p10[1]-p01[1], p10[2]-p01[2]);
    Vec3 raw = d1.cross(d2);
    return raw.normalized();
  };

  Vec3 base_normal = compute_normal(n00, n10, n01, n11);

  // Check each corner node, z-direction
  for (int corner = 0; corner < 4; ++corner) {
    double* node = nullptr;
    double* d_node = nullptr;
    switch (corner) {
      case 0: node = n00; d_node = d_n00; break;
      case 1: node = n10; d_node = d_n10; break;
      case 2: node = n01; d_node = d_n01; break;
      case 3: node = n11; d_node = d_n11; break;
    }

    for (int dir = 0; dir < 3; ++dir) {
      double saved = node[dir];

      node[dir] = saved + h;
      Vec3 n_plus = compute_normal(n00, n10, n01, n11);

      node[dir] = saved - h;
      Vec3 n_minus = compute_normal(n00, n10, n01, n11);

      node[dir] = saved;

      // FD of (d_normal . normal) w.r.t. node[dir]
      double fd = (d_normal[0] * (n_plus[0] - n_minus[0]) +
                   d_normal[1] * (n_plus[1] - n_minus[1]) +
                   d_normal[2] * (n_plus[2] - n_minus[2])) / (2.0 * h);

      EXPECT_NEAR(d_node[dir], fd, 1e-5)
          << "Normal gradient mismatch at corner=" << corner
          << " dir=" << dir;
    }
  }
}

}  // namespace
}  // namespace ad
}  // namespace vlm
