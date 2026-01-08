// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_VLM_H_
#define VLM_VLM_H_

/// @file vlm.h
/// @brief Main facade header for the VLM (Vortex Lattice Method) library.
///
/// This header includes all public API components needed for typical usage.
/// For specialized needs, include individual headers directly.
///
/// @section quick_start Quick Start
///
/// @code
/// #include <vlm/vlm.h>
///
/// int main() {
///   // Load aircraft geometry from AVL file
///   vlm::AvlReader reader;
///   auto aircraft = reader.Read("aircraft.avl");
///
///   // Create solver and set up geometry
///   vlm::VortexLatticeSolver solver;
///   solver.SetAircraft(aircraft);
///
///   // Solve at a flight condition
///   vlm::FlightCondition condition;
///   condition.alpha = 5.0 * vlm::kDegreesToRadians;
///   condition.mach = 0.0;
///
///   solver.Solve(condition);
///
///   // Compute forces
///   vlm::ForceCalculator calculator;
///   auto results = calculator.Compute(solver, condition);
///
///   std::cout << "CL = " << results.CL << std::endl;
///   std::cout << "CD = " << results.CD << std::endl;
/// }
/// @endcode

// Core types and constants
#include "vlm/core/types.h"

// Geometry
#include "vlm/geometry/aircraft.h"
#include "vlm/geometry/airfoil.h"
#include "vlm/geometry/body.h"
#include "vlm/geometry/section.h"
#include "vlm/geometry/surface.h"

// I/O
#include "vlm/io/avl_reader.h"

// Solver
#include "vlm/solver/horseshoe_vortex.h"
#include "vlm/solver/influence_matrix.h"
#include "vlm/solver/vortex_lattice.h"
#include "vlm/solver/vortex_lattice_solver.h"

// Aerodynamics
#include "vlm/aero/aerodynamic_results.h"
#include "vlm/aero/force_calculator.h"
#include "vlm/aero/stability_derivatives.h"

// Trim
#include "vlm/trim/flight_constraints.h"
#include "vlm/trim/trim_solver.h"

// Dynamics
#include "vlm/dynamics/eigenmode_analysis.h"
#include "vlm/dynamics/mass_properties.h"

// Utilities
#include "vlm/utils/cubic_spline.h"
#include "vlm/utils/math_utilities.h"

/// @namespace vlm
/// @brief Vortex Lattice Method library namespace.
///
/// The VLM library provides a complete implementation of the Vortex Lattice
/// Method for aerodynamic analysis of aircraft configurations. Key features:
///
/// - **Geometry**: Support for lifting surfaces, bodies, and control surfaces
/// - **AVL Compatibility**: Read AVL geometry files directly
/// - **Solver**: Efficient horseshoe vortex implementation with Prandtl-Glauert
/// - **Forces**: Near-field Kutta-Joukowski and far-field Trefftz plane methods
/// - **Stability**: Full stability derivative computation
/// - **Trim**: Newton-Raphson trim solver with flexible constraints
/// - **Dynamics**: Eigenmode analysis for dynamic stability
/// - **AD**: Automatic differentiation for gradient-based optimization

#endif  // VLM_VLM_H_
