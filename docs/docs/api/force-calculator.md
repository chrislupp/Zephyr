---
sidebar_position: 6
---

# ForceCalculator

```cpp
#include "pegasus/aero/force_calculator.h"
```

Computes aerodynamic forces and moments from a VLM solution using the Kutta-Joukowski theorem and Trefftz plane analysis.

## Methods

```cpp
AerodynamicResults ComputeForces(
    const Aircraft& aircraft,
    const VortexLattice& mesh,
    const SolutionState& solution,
    const FlightCondition& condition);
```

Computes all aerodynamic forces and moments. Uses near-field Kutta-Joukowski for lift and side force, and Trefftz plane for induced drag.

### Options

```cpp
void SetIncludeInducedVelocity(bool include);
void SetTrefftzCoreRadius(Scalar chord_fraction, Scalar width_fraction);
```

## AerodynamicResults

```cpp
struct AerodynamicResults {
    // Force coefficients
    Scalar CL;    // lift coefficient
    Scalar CD;    // total drag coefficient
    Scalar CY;    // side force coefficient
    Scalar CDi;   // induced drag coefficient (Trefftz plane)

    // Moment coefficients
    Scalar Cl;    // rolling moment coefficient
    Scalar Cm;    // pitching moment coefficient
    Scalar Cn;    // yawing moment coefficient

    // Per-surface breakdowns
    // ...
};
```
