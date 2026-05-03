---
sidebar_position: 11
---

# GradientCalculator

```cpp
#include "vlm/ad/gradient_calculator.h"
```

High-level interface for computing analytical gradients of aerodynamic coefficients using the adjoint method.

## Setup

```cpp
void SetAircraft(const Aircraft& aircraft);
void SetReferencePoint(const Vec3& ref);
void SetMesh(const VortexLattice& mesh);
void SetFiniteDiffStep(Scalar step);
```

## Flight Condition Derivatives

```cpp
FullGradients ComputeFlightConditionDerivatives();
AeroGradients ComputeAlphaDerivatives();
AeroGradients ComputeBetaDerivatives();
AeroGradients ComputeRollRateDerivatives();
AeroGradients ComputePitchRateDerivatives();
AeroGradients ComputeYawRateDerivatives();
AeroGradients ComputeMachDerivatives();
```

### AeroGradients

```cpp
struct AeroGradients {
    Scalar dCL, dCD, dCY;     // force coefficient derivatives
    Scalar dCl, dCm, dCn;     // moment coefficient derivatives
    Scalar dCDi;               // induced drag derivative
};
```

### FullGradients

```cpp
struct FullGradients {
    AeroGradients alpha_derivs;
    AeroGradients beta_derivs;
    AeroGradients p_derivs;
    AeroGradients q_derivs;
    AeroGradients r_derivs;
    AeroGradients mach_derivs;
};
```

## Geometry Gradients

```cpp
GeometryGradients ComputeGeometryGradients();
SectionGeometryGradients ComputeSectionGradient(
    Index surface_idx, Index section_idx, GeometryVariable variable);
```

### GeometryVariable

```cpp
enum class GeometryVariable {
    kChord, kTwist, kXle, kYle, kZle
};
```

### SectionGeometryGradients

```cpp
struct SectionGeometryGradients {
    std::string surface_name;
    Index section_index;
    AeroGradients chord_derivs;
    AeroGradients twist_derivs;
    AeroGradients xle_derivs;
    AeroGradients yle_derivs;
    AeroGradients zle_derivs;
};
```

## Mesh Node Gradients

For direct mesh-level sensitivity (useful in aerostructural coupling):

```cpp
auto mesh_grads = ComputeMeshNodeGradients();
```

Returns derivatives of aerodynamic coefficients with respect to each mesh node's 3D position.
