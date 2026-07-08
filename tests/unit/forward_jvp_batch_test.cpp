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

/// @file forward_jvp_batch_test.cpp
/// @brief Batched forward JVPs must match the per-direction functions.

#include "pegasus/aero/force_calculator.h"
#include "pegasus/gradients/forward_jvp.h"
#include "pegasus/solver/vortex_lattice_solver.h"

#include <random>

#include <gtest/gtest.h>

namespace pegasus {
namespace gradients {
namespace {

/// Mildly twisted rectangular wing so all derivative paths are exercised
SurfaceMesh BuildTwistedMesh(Index nx, Index ny, Scalar chord, Scalar span) {
    SurfaceMesh sm;
    sm.name = "Wing";
    sm.nx = nx;
    sm.ny = ny;
    sm.nodes.resize((nx + 1) * (ny + 1), 3);

    for (Index i = 0; i <= nx; ++i) {
        for (Index j = 0; j <= ny; ++j) {
            Scalar x = chord * Scalar(i) / Scalar(nx);
            Scalar y = span * Scalar(j) / Scalar(ny);
            Scalar theta = 0.05 * y / span;
            Scalar xr = x - 0.25 * chord;
            sm.SetNode(i, j,
                       Vec3(0.25 * chord + xr * std::cos(theta), y,
                            -xr * std::sin(theta)));
        }
    }
    return sm;
}

TEST(ForwardJvpBatchTest, BatchMatchesPerDirectionJVPs) {
    const Index nx = 3;
    const Index ny = 6;
    const Index n = nx * ny;
    const Index kk = 7;

    SurfaceMesh sm = BuildTwistedMesh(nx, ny, 1.0, 4.0);
    MeshInput input;
    input.surfaces.push_back(sm);
    input.reference.reference_area = 4.0;
    input.reference.reference_chord = 1.0;
    input.reference.reference_span = 4.0;

    FlightCondition condition;
    condition.alpha = 0.08;

    VortexLatticeSolver solver;
    solver.SetMesh(input);
    SolutionState solution = solver.Solve(condition);
    ASSERT_TRUE(solution.is_valid);

    // Random node perturbation directions
    std::mt19937 gen(1234);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);

    std::vector<VortexGeometryPerturbation> perts(kk);
    VortexGeometryPerturbationBatch batch;
    batch.Resize(n, kk);
    for (Index k = 0; k < kk; ++k) {
        MeshNodeGradients delta_nodes(nx, ny);
        for (Index d = 0; d < delta_nodes.d_nodes.size(); ++d) {
            delta_nodes.d_nodes[d] = dist(gen);
        }
        perts[k] =
            PropagateNodePerturbationsToVortex(sm, 0, nx, ny, delta_nodes);
        batch.SetColumn(k, perts[k]);
    }

    // Residual JVPs
    MatX delta_R_batch =
        ComputeResidualGeometryJVPBatch(solver, solution, condition, batch);
    ASSERT_EQ(delta_R_batch.rows(), n);
    ASSERT_EQ(delta_R_batch.cols(), kk);

    MatX delta_gammas(n, kk);
    for (Index k = 0; k < kk; ++k) {
        VecX delta_R =
            ComputeResidualGeometryJVP(solver, solution, condition, perts[k]);
        EXPECT_LT((delta_R_batch.col(k) - delta_R).norm(),
                  1e-12 * std::max(1.0, delta_R.norm()))
            << "residual JVP mismatch in direction " << k;
        delta_gammas.col(k) = -solver.SolveForward(delta_R);
    }

    // Matrix SolveForward overload matches column solves
    MatX delta_gammas_mat = -solver.SolveForward(delta_R_batch);
    EXPECT_LT((delta_gammas_mat - delta_gammas).norm(),
              1e-12 * std::max(1.0, delta_gammas.norm()));

    // Force JVPs
    MatX delta_F_batch = ComputePerVortexForceTotalJVPBatch(
        solver, solution, condition, batch, delta_gammas);
    ASSERT_EQ(delta_F_batch.rows(), 3 * n);
    ASSERT_EQ(delta_F_batch.cols(), kk);

    for (Index k = 0; k < kk; ++k) {
        VecX delta_F = ComputePerVortexForceTotalJVP(
            solver, solution, condition, perts[k], delta_gammas.col(k));
        EXPECT_LT((delta_F_batch.col(k) - delta_F).norm(),
                  1e-12 * std::max(1.0, delta_F.norm()))
            << "force JVP mismatch in direction " << k;
    }
}

/// Inner product between a forward geometry perturbation and a reverse
/// geometry adjoint seed, summed over all four geometry quantities.
double GeomDot(const VortexGeometryPerturbation& v, const MeshGradients& g) {
    double s = 0.0;
    for (Eigen::Index k = 0; k < v.d_rv1.size(); ++k) {
        s += v.d_rv1(k) * g.d_rv1[size_t(k)];
        s += v.d_rv2(k) * g.d_rv2[size_t(k)];
        s += v.d_rc(k) * g.d_rc[size_t(k)];
        s += v.d_normal(k) * g.d_normal[size_t(k)];
    }
    return s;
}

TEST(ForwardJvpBatchTest, GeometryVJPsAreExactTransposes) {
    const Index nx = 3;
    const Index ny = 6;
    const Index n = nx * ny;

    SurfaceMesh sm = BuildTwistedMesh(nx, ny, 1.0, 4.0);
    MeshInput input;
    input.surfaces.push_back(sm);
    input.reference.reference_area = 4.0;
    input.reference.reference_chord = 1.0;
    input.reference.reference_span = 4.0;

    FlightCondition condition;
    condition.alpha = 0.08;

    VortexLatticeSolver solver;
    solver.SetMesh(input);
    SolutionState solution = solver.Solve(condition);
    ASSERT_TRUE(solution.is_valid);

    std::mt19937 gen(99);
    std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);

    // Random forward geometry perturbation direction v.
    MeshNodeGradients delta_nodes(nx, ny);
    for (Index d = 0; d < delta_nodes.d_nodes.size(); ++d) {
        delta_nodes.d_nodes[d] = dist(gen);
    }
    VortexGeometryPerturbation v =
        PropagateNodePerturbationsToVortex(sm, 0, nx, ny, delta_nodes);

    // --- Residual: <bar_R, J_R v> == <VJP_R(bar_R), v> ---
    VecX bar_R(n);
    for (Index i = 0; i < n; ++i)
        bar_R(i) = dist(gen);
    VecX dR = ComputeResidualGeometryJVP(solver, solution, condition, v);
    MeshGradients gR = ComputeResidualGeometryVJP(solver, solution, condition,
                                                  bar_R);
    double lhs_R = bar_R.dot(dR);
    double rhs_R = GeomDot(v, gR);
    EXPECT_NEAR(lhs_R, rhs_R, 1e-9 * std::max(1.0, std::abs(lhs_R)));

    // --- Force geometry: <bar_F, J_F,geom v> == <VJP_F,geom(bar_F), v> ---
    // Isolate the geometry terms of the force JVP with delta_gamma = 0.
    VecX zero_gamma = VecX::Zero(n);
    VecX dF = ComputePerVortexForceTotalJVP(solver, solution, condition, v,
                                            zero_gamma);
    VecX bar_F(3 * n);
    for (Index i = 0; i < 3 * n; ++i)
        bar_F(i) = dist(gen);
    MeshGradients gF = ComputePerVortexForceGeometryVJP(solver, solution,
                                                        condition, bar_F);
    double lhs_F = bar_F.dot(dF);
    double rhs_F = GeomDot(v, gF);
    EXPECT_NEAR(lhs_F, rhs_F, 1e-9 * std::max(1.0, std::abs(lhs_F)));
}

}  // namespace
}  // namespace gradients
}  // namespace pegasus
