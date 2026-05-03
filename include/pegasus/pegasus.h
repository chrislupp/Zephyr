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

#ifndef PEGASUS_PEGASUS_H_
#define PEGASUS_PEGASUS_H_

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
///   pegasus::AvlReader reader;
///   auto aircraft = reader.Read("aircraft.avl");
///
///   // Create solver and set up geometry
///   pegasus::VortexLatticeSolver solver;
///   solver.SetAircraft(aircraft);
///
///   // Solve at a flight condition
///   pegasus::FlightCondition condition;
///   condition.alpha = 5.0 * pegasus::kDegreesToRadians;
///   condition.mach = 0.0;
///
///   solver.Solve(condition);
///
///   // Compute forces
///   pegasus::ForceCalculator calculator;
///   auto results = calculator.Compute(solver, condition);
///
///   std::cout << "CL = " << results.CL << std::endl;
///   std::cout << "CD = " << results.CD << std::endl;
/// }
/// @endcode

// Core types and constants
#include "pegasus/core/types.h"

// Geometry
#include "pegasus/geometry/aircraft.h"
#include "pegasus/geometry/airfoil.h"
#include "pegasus/geometry/body.h"
#include "pegasus/geometry/section.h"
#include "pegasus/geometry/surface.h"

// I/O
#include "pegasus/io/avl_reader.h"

// Solver
#include "pegasus/solver/horseshoe_vortex.h"
#include "pegasus/solver/influence_matrix.h"
#include "pegasus/solver/vortex_lattice.h"
#include "pegasus/solver/vortex_lattice_solver.h"

// Aerodynamics
#include "pegasus/aero/aerodynamic_results.h"
#include "pegasus/aero/force_calculator.h"
#include "pegasus/aero/stability_derivatives.h"

// Trim
#include "pegasus/trim/flight_constraints.h"
#include "pegasus/trim/trim_solver.h"

// Dynamics
#include "pegasus/dynamics/eigenmode_analysis.h"
#include "pegasus/dynamics/mass_properties.h"

// Utilities
#include "pegasus/utils/cubic_spline.h"
#include "pegasus/utils/math_utilities.h"

// Visualization (VTK export - optional, check VtkExporter::IsAvailable())
#include "pegasus/visualization/vtk_export.h"

/// @namespace pegasus
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

#endif  // PEGASUS_PEGASUS_H_
