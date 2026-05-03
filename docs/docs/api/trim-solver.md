---
sidebar_position: 8
---

# TrimSolver

```cpp
#include "pegasus/trim/trim_solver.h"
```

Newton-Raphson trim solver for finding aircraft equilibrium at specified flight constraints.

## Methods

```cpp
void SetOptions(const TrimOptions& options);
TrimOptions GetOptions() const;

TrimResult Solve(
    const Aircraft& aircraft,
    VortexLatticeSolver& solver,
    const TrimConstraints& constraints,
    const OperatingVariables& initial_state = {});
```

## TrimOptions

```cpp
struct TrimOptions {
    int max_iterations = 50;
    Scalar tolerance = 1e-8;
    Scalar max_alpha = 30.0 * kDegreesToRadians;
    Scalar max_beta = 30.0 * kDegreesToRadians;
    Scalar max_rate = 1.0;
    Scalar relaxation = 1.0;
    bool verbose = false;
};
```

## TrimResult

```cpp
struct TrimResult {
    bool converged;
    int iterations;
    Scalar residual_norm;
    OperatingVariables state;
    AerodynamicResults aero;
    StabilityDerivatives derivatives;
    std::string error_message;
};
```

## Presets

```cpp
TrimConstraints TrimPresets::LevelFlight(Scalar target_CL);
```

Creates constraints for level, unaccelerated flight at a specified $C_L$.
