---
sidebar_position: 2
---

# Quick Start

This guide walks through a minimal VLM analysis: defining a rectangular wing, solving for the circulation distribution, and computing aerodynamic forces.

## Include the Library

VLM provides a single facade header that includes all public API components:

```cpp
#include "pegasus/vlm.h"

using namespace pegasus;
```

## Define the Geometry

Create a rectangular wing with a semi-span of 5.0 and chord of 1.0:

```cpp
// Define a lifting surface
Surface wing("Wing");
wing.SetChordwisePanelCount(6);
wing.SetSpanwisePanelCount(20);
wing.SetYDuplicate(0.0);  // Mirror about Y=0

// Add root and tip sections
Section root(Vec3(0.0, 0.0, 0.0), 1.0);  // chord = 1.0
Section tip(Vec3(0.0, 5.0, 0.0), 1.0);   // semi-span = 5.0
wing.AddSection(root);
wing.AddSection(tip);

// Create the aircraft configuration
Aircraft aircraft("Rectangular Wing");
aircraft.AddSurface(wing);
aircraft.SetReferenceArea(10.0);
aircraft.SetReferenceChord(1.0);
aircraft.SetReferenceSpan(10.0);
```

## Solve

Set up the solver, define the flight condition, and solve:

```cpp
VortexLatticeSolver solver;
solver.SetAircraft(aircraft);

FlightCondition condition;
condition.alpha = 5.0 * kDegreesToRadians;

SolutionState solution = solver.Solve(condition);
```

## Compute Forces

```cpp
ForceCalculator calculator;
AerodynamicResults results = calculator.ComputeForces(
    aircraft, solver.GetMesh(), solution, condition);

std::cout << "CL = " << results.CL << std::endl;
std::cout << "CD = " << results.CD << std::endl;
std::cout << "L/D = " << results.CL / results.CD << std::endl;
```

## Loading AVL Files

VLM can parse Athena Vortex Lattice geometry files:

```cpp
AvlReader reader;
Aircraft aircraft = reader.ReadFile("aircraft.avl");

VortexLatticeSolver solver;
solver.SetAircraft(aircraft);

FlightCondition condition;
condition.alpha = 5.0 * kDegreesToRadians;

SolutionState solution = solver.Solve(condition);

ForceCalculator calculator;
AerodynamicResults results = calculator.ComputeForces(
    aircraft, solver.GetMesh(), solution, condition);

std::cout << "Aircraft: " << aircraft.GetTitle() << std::endl;
std::cout << "CL = " << results.CL << std::endl;
```

## Next Steps

- Read the [Solver Workflow](../guide/solver-workflow.md) guide for a deeper understanding of the analysis pipeline
- Explore [Stability Analysis](../guide/stability-analysis.md) for computing derivatives and dynamic modes
- Check the [Examples](../examples/basic-analysis.md) for complete programs
