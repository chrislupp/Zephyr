---
sidebar_position: 9
---

# EigenmodeAnalysis

```cpp
#include "pegasus/dynamics/eigenmode_analysis.h"
```

Linearized dynamic stability analysis through eigenvalue decomposition of the equations of motion.

## Methods

```cpp
EigenmodeResults Analyze(
    const Aircraft& aircraft,
    const MassProperties& mass,
    const StabilityDerivatives& derivs,
    const FlightState& flight);

EigenmodeResults AnalyzeLongitudinal(
    const Aircraft& aircraft,
    const MassProperties& mass,
    const StabilityDerivatives& derivs,
    const FlightState& flight);

EigenmodeResults AnalyzeLateralDirectional(
    const Aircraft& aircraft,
    const MassProperties& mass,
    const StabilityDerivatives& derivs,
    const FlightState& flight);
```

## MassProperties

```cpp
struct MassProperties {
    Scalar mass;         // total mass (kg)
    Scalar Ixx, Iyy, Izz;  // moments of inertia (kg*m^2)
    Scalar Ixz;          // product of inertia
};
```

## FlightState

```cpp
struct FlightState {
    Scalar velocity;   // airspeed (m/s)
    Scalar density;    // air density (kg/m^3)
};
```

## Eigenmode

```cpp
struct Eigenmode {
    std::complex<Scalar> eigenvalue;
    ModeType type;

    Scalar GetDamping() const;
    Scalar GetFrequency() const;
    Scalar GetPeriod() const;
    Scalar GetTimeConstant() const;
    Scalar GetDampingRatio() const;
    bool IsStable() const;
    bool IsOscillatory() const;
};
```

## ModeType

```cpp
enum class ModeType {
    kShortPeriod,
    kPhugoid,
    kDutchRoll,
    kRollSubsidence,
    kSpiralMode,
    kOther
};
```

## EigenmodeResults

```cpp
struct EigenmodeResults {
    bool is_valid;
    std::vector<Eigenmode> modes;
    MatX system_matrix;

    std::vector<const Eigenmode*> FindModesByType(ModeType type) const;
    const Eigenmode* GetShortPeriod() const;
    const Eigenmode* GetPhugoid() const;
    const Eigenmode* GetDutchRoll() const;
    bool IsStable() const;
};
```
