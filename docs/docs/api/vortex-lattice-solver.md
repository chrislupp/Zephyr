---
sidebar_position: 5
---

# VortexLatticeSolver

```cpp
#include "pegasus/solver/vortex_lattice_solver.h"
```

Core VLM solver that builds the aerodynamic influence coefficient (AIC) matrix and solves for circulation distributions.

## Setup

```cpp
void SetAircraft(const Aircraft& aircraft);
void SetMesh(const VortexLattice& mesh);
void SetSymmetry(SymmetryType symmetry);
void SetReferencePoint(const Vec3& point);
```

## Solving

```cpp
SolutionState Solve(const FlightCondition& condition);
```

Solves the linear system $A \Gamma = b$ for the circulation distribution. The AIC matrix is built lazily on the first call.

```cpp
void RebuildAIC();
void SolveUnitResponses();
bool IsAICValid() const;
bool AreUnitResponsesValid() const;
```

## Accessing Results

```cpp
const SolutionState& GetSolution() const;
const VortexLattice& GetMesh() const;
const MatX& GetAICMatrix() const;
```

## Key Types

### FlightCondition

```cpp
struct FlightCondition {
    Scalar alpha = 0.0;       // angle of attack (rad)
    Scalar beta = 0.0;        // sideslip angle (rad)
    Scalar roll_rate = 0.0;   // non-dimensional roll rate (pb/2V)
    Scalar pitch_rate = 0.0;  // non-dimensional pitch rate (qc/2V)
    Scalar yaw_rate = 0.0;    // non-dimensional yaw rate (rb/2V)
    Scalar mach = 0.0;        // Mach number
    Scalar dynamic_pressure = 0.5;
    std::map<std::string, Scalar> control_deflections;
};
```

### SolutionState

```cpp
struct SolutionState {
    VecX circulation;             // circulation per panel
    VecX circulation_per_u;       // unit response per velocity component
    VecX circulation_per_v;
    VecX circulation_per_w;
    VecX circulation_per_p;
    VecX circulation_per_q;
    VecX circulation_per_r;
    VecX induced_velocity;
    bool is_valid = false;
};
```
