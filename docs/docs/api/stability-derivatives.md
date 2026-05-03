---
sidebar_position: 7
---

# StabilityDerivativeCalculator

```cpp
#include "vlm/aero/stability_derivatives.h"
```

Computes stability and control derivatives by evaluating unit perturbation responses.

## Methods

```cpp
StabilityDerivatives ComputeDerivatives(
    const Aircraft& aircraft,
    const VortexLatticeSolver& solver,
    const FlightCondition& condition);

BodyAxisDerivatives ComputeBodyAxisDerivatives(
    const Aircraft& aircraft,
    const VortexLatticeSolver& solver,
    const FlightCondition& condition);
```

## StabilityDerivatives

Contains derivatives organized by perturbation variable:

### Angle of Attack Derivatives
`CL_alpha`, `CD_alpha`, `CY_alpha`, `Cl_alpha`, `Cm_alpha`, `Cn_alpha`

### Sideslip Derivatives
`CL_beta`, `CD_beta`, `CY_beta`, `Cl_beta`, `Cm_beta`, `Cn_beta`

### Rate Derivatives
- Roll rate: `CL_p`, `Cl_p`, `Cn_p`, ...
- Pitch rate: `CL_q`, `Cm_q`, ...
- Yaw rate: `CY_r`, `Cl_r`, `Cn_r`, ...

### Control Derivatives

```cpp
struct ControlDerivatives {
    std::string name;
    Scalar deflection;
    Scalar gain;
    Scalar CL_d, CY_d, CD_d;
    Scalar Cl_d, Cm_d, Cn_d;
};
```

Accessed via `derivs.control_derivatives` (one per control surface).

### Static Stability

```cpp
Scalar neutral_point_x;  // neutral point x-coordinate
Scalar static_margin;     // static margin (fraction of reference chord)
```
