---
sidebar_position: 3
---

# Stability Modes

This example demonstrates full stability derivative computation and eigenmode analysis for a three-surface aircraft.

## Overview

The program:
1. Creates a wing-tail-fin aircraft with sweep and dihedral
2. Computes longitudinal and lateral-directional stability derivatives
3. Performs eigenmode analysis to identify dynamic modes
4. Reports mode characteristics: frequency, damping ratio, and stability

## Code

```cpp
#include "pegasus/vlm.h"

#include <iostream>

using namespace pegasus;

int main() {
    // Build a three-surface aircraft (wing, htail, vtail)
    // ... (geometry setup as in stability_analysis example) ...

    VortexLatticeSolver solver;
    solver.SetAircraft(aircraft);

    FlightCondition condition;
    condition.alpha = 5.0 * kDegreesToRadians;
    solver.Solve(condition);

    // Compute stability derivatives
    StabilityDerivativeCalculator deriv_calc;
    StabilityDerivatives derivs = deriv_calc.ComputeDerivatives(
        aircraft, solver, condition);

    std::cout << "CL_alpha = " << derivs.CL_alpha << " /rad" << std::endl;
    std::cout << "Cm_alpha = " << derivs.Cm_alpha << " /rad" << std::endl;
    std::cout << "Static margin = " << derivs.static_margin << std::endl;

    // Mass properties
    MassProperties mass;
    mass.mass = 1000.0;
    mass.Ixx = 1000.0;
    mass.Iyy = 2000.0;
    mass.Izz = 2500.0;
    mass.Ixz = 100.0;

    FlightState flight;
    flight.velocity = 50.0;
    flight.density = 1.225;

    // Longitudinal eigenmode analysis
    EigenmodeAnalysis analysis;
    EigenmodeResults lon = analysis.AnalyzeLongitudinal(
        aircraft, mass, derivs, flight);

    std::cout << "\nLongitudinal modes:" << std::endl;
    for (const auto& mode : lon.modes) {
        std::cout << "  " << GetModeTypeName(mode.type) << ": ";
        if (mode.IsOscillatory()) {
            std::cout << "freq=" << mode.GetFrequency() << " rad/s"
                      << " damping=" << mode.GetDampingRatio();
        } else {
            std::cout << "time_const=" << mode.GetTimeConstant() << " s";
        }
        std::cout << (mode.IsStable() ? " [STABLE]" : " [UNSTABLE]")
                  << std::endl;
    }

    // Lateral-directional eigenmode analysis
    EigenmodeResults lat = analysis.AnalyzeLateralDirectional(
        aircraft, mass, derivs, flight);

    std::cout << "\nLateral-directional modes:" << std::endl;
    for (const auto& mode : lat.modes) {
        std::cout << "  " << GetModeTypeName(mode.type) << ": ";
        if (mode.IsOscillatory()) {
            std::cout << "freq=" << mode.GetFrequency() << " rad/s"
                      << " damping=" << mode.GetDampingRatio();
        } else {
            std::cout << "time_const=" << mode.GetTimeConstant() << " s";
        }
        std::cout << (mode.IsStable() ? " [STABLE]" : " [UNSTABLE]")
                  << std::endl;
    }

    return 0;
}
```

## Running

```bash
./build/examples/stability_analysis
```

## Dynamic Modes

### Longitudinal
- **Short Period** — fast pitch oscillation (typically 1-5 rad/s)
- **Phugoid** — slow speed/altitude exchange (typically 0.1-0.5 rad/s)

### Lateral-Directional
- **Dutch Roll** — coupled yaw-roll oscillation
- **Roll Subsidence** — pure roll damping (non-oscillatory)
- **Spiral Mode** — slow heading divergence/convergence (non-oscillatory)
