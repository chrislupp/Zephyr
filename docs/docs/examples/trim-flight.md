---
sidebar_position: 2
---

# Trim Flight

This example demonstrates the Newton-Raphson trim solver for finding aircraft equilibrium in level flight.

## Overview

The program:
1. Creates a two-surface aircraft (wing + horizontal tail)
2. Trims the aircraft at multiple target $C_L$ values (0.3, 0.5, 0.7, 1.0)
3. Reports the required angle of attack, resulting $C_m$, and convergence iterations

## Code

```cpp
#include "vlm/vlm.h"

#include <iostream>

using namespace vlm;

int main() {
    // Main wing
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(20);
    wing.SetYDuplicate(0.0);

    Section wing_root(Vec3(0.0, 0.0, 0.0), 2.0);
    Section wing_tip(Vec3(0.5, 7.5, 0.0), 1.0);
    wing.AddSection(wing_root);
    wing.AddSection(wing_tip);

    // Horizontal tail
    Surface tail("Tail");
    tail.SetChordwisePanelCount(4);
    tail.SetSpanwisePanelCount(10);
    tail.SetYDuplicate(0.0);

    Section tail_root(Vec3(5.0, 0.0, 0.0), 1.0);
    Section tail_tip(Vec3(5.2, 2.5, 0.0), 0.7);
    tail.AddSection(tail_root);
    tail.AddSection(tail_tip);

    Aircraft aircraft("Trim Example");
    aircraft.AddSurface(wing);
    aircraft.AddSurface(tail);
    aircraft.SetReferenceArea(30.0);
    aircraft.SetReferenceChord(2.0);
    aircraft.SetReferenceSpan(15.0);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    TrimSolver trim_solver;

    // Trim at multiple CL targets
    for (double cl_target : {0.3, 0.5, 0.7, 1.0}) {
        TrimConstraints constraints = TrimPresets::LevelFlight(cl_target);
        TrimResult result = trim_solver.Solve(aircraft, solver, constraints);

        if (result.converged) {
            std::cout << "CL_target=" << cl_target
                      << "  alpha=" << result.state.alpha * kRadiansToDegrees << " deg"
                      << "  CL=" << result.aero.CL
                      << "  Cm=" << result.aero.Cm
                      << "  iters=" << result.iterations << std::endl;
        }
    }

    return 0;
}
```

## Running

```bash
./build/examples/trim_flight
```
