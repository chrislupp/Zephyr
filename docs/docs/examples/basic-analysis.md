---
sidebar_position: 1
---

# Basic Analysis

This example demonstrates a basic aerodynamic analysis of a rectangular wing, sweeping through angles of attack and computing lift, drag, and efficiency.

## Overview

The program:
1. Creates a rectangular wing with specified panel density
2. Sweeps angle of attack from $-2°$ to $10°$
3. Computes $C_L$, $C_D$, and $L/D$ at each angle
4. Compares the computed lift curve slope with lifting line theory

## Code

```cpp
#include "pegasus/vlm.h"

#include <iomanip>
#include <iostream>

using namespace pegasus;

int main() {
    // Create a rectangular wing
    Surface wing("Wing");
    wing.SetChordwisePanelCount(6);
    wing.SetSpanwisePanelCount(20);
    wing.SetYDuplicate(0.0);

    Section root(Vec3(0.0, 0.0, 0.0), 1.0);
    Section tip(Vec3(0.0, 5.0, 0.0), 1.0);
    wing.AddSection(root);
    wing.AddSection(tip);

    Aircraft aircraft("Rectangular Wing");
    aircraft.AddSurface(wing);
    aircraft.SetReferenceArea(10.0);
    aircraft.SetReferenceChord(1.0);
    aircraft.SetReferenceSpan(10.0);

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    ForceCalculator calculator;

    // Angle of attack sweep
    std::cout << std::setw(10) << "Alpha"
              << std::setw(12) << "CL"
              << std::setw(12) << "CD"
              << std::setw(12) << "L/D" << std::endl;

    for (double alpha_deg = -2.0; alpha_deg <= 10.0; alpha_deg += 1.0) {
        FlightCondition condition;
        condition.alpha = alpha_deg * kDegreesToRadians;

        SolutionState solution = solver.Solve(condition);
        AerodynamicResults results = calculator.ComputeForces(
            aircraft, solver.GetMesh(), solution, condition);

        std::cout << std::setw(10) << alpha_deg
                  << std::setw(12) << results.CL
                  << std::setw(12) << results.CD
                  << std::setw(12) << results.CL / results.CD
                  << std::endl;
    }

    return 0;
}
```

## Running

```bash
./build/examples/basic_analysis
```

## Expected Output

The lift curve slope should be close to the theoretical lifting line value for an aspect ratio 10 wing:

$$
C_{L_\alpha} = \frac{2\pi}{1 + 2/AR} \approx 5.24 \text{ /rad}
$$
