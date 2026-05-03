// Single-analysis scaling benchmark: native C++ VLM (no Philote/gRPC).
//
// Sweeps spanwise panel count with fixed nx=6, runs mesh generation + solve
// multiple times per data point for timing accuracy. Outputs CSV to stdout.
//
// Build:  make benchmark_analysis
// Run:    ./benchmark_analysis > results_vlm.csv

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include "vlm/vlm.h"

int main() {
  constexpr double span = 10.0;
  constexpr double chord = 1.0;
  constexpr double alpha_rad = 5.0 * M_PI / 180.0;
  constexpr int nx = 6;
  constexpr int n_repeats = 50;

  std::vector<int> ny_values = {8, 12, 16, 24, 32, 48, 64, 96};

  std::cout << "ny,panels,time_ms" << std::endl;

  for (int ny : ny_values) {
    int n_stations = ny + 1;

    // Build mesh input (matches the Philote benchmark geometry)
    vlm::MeshInput mesh_input;
    mesh_input.reference.reference_area = span * chord;
    mesh_input.reference.reference_chord = chord;
    mesh_input.reference.reference_span = span;
    mesh_input.reference.moment_reference_point = {0.25 * chord, 0.0, 0.0};

    vlm::SurfaceMesh surface;
    surface.name = "wing";
    surface.nx = nx;
    surface.ny = ny;
    surface.has_wake = true;
    surface.y_duplicated = true;
    surface.y_duplicate_plane = 0.0;
    surface.affected_by_freestream = true;
    surface.nodes.resize(surface.NumNodes(), 3);

    // Generate mesh nodes: cosine chordwise, linear spanwise
    for (int i = 0; i <= nx; ++i) {
      double cf = 0.5 * (1.0 - std::cos(M_PI * i / nx));
      for (int j = 0; j < n_stations; ++j) {
        double y = (static_cast<double>(j) / ny) * (span / 2.0);
        surface.SetNode(i, j, {cf * chord, y, 0.0});
      }
    }
    mesh_input.surfaces.push_back(surface);

    vlm::FlightCondition condition;
    condition.alpha = alpha_rad;
    condition.dynamic_pressure = 0.5 * 1.225 * 25.0 * 25.0;

    // Warm-up run
    {
      vlm::VortexLatticeSolver solver;
      solver.SetMesh(mesh_input);
      solver.RebuildAIC();
      solver.Solve(condition);
    }

    // Timed runs
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int r = 0; r < n_repeats; ++r) {
      vlm::VortexLatticeSolver solver;
      solver.SetMesh(mesh_input);
      solver.RebuildAIC();
      vlm::SolutionState state = solver.Solve(condition);

      vlm::ForceCalculator calc;
      calc.ComputeForces(mesh_input.reference, solver.GetMesh(),
                         state, condition);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    double per_run_ms = elapsed_ms / n_repeats;

    int total_panels = nx * 2 * ny;  // full-span
    std::cout << ny << "," << total_panels << "," << per_run_ms << std::endl;
    std::cerr << "ny=" << ny << " (" << total_panels << " panels): "
              << per_run_ms << " ms/run" << std::endl;
  }

  return 0;
}
