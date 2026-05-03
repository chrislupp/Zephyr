---
sidebar_position: 4
---

# Section

```cpp
#include "vlm/geometry/section.h"
```

Defines a spanwise station on a lifting surface, specifying chord, twist, airfoil, and control surfaces.

## Constructors

```cpp
Section();
Section(const Vec3& leading_edge_position,
        Scalar chord,
        Scalar twist_angle = 0.0,
        std::shared_ptr<Airfoil> airfoil = nullptr);
```

**Parameters:**
- `leading_edge_position` — 3D position of the section leading edge
- `chord` — chord length
- `twist_angle` — geometric twist in degrees (positive = leading edge up)
- `airfoil` — optional airfoil definition

## Geometry

```cpp
void SetLeadingEdgePosition(const Vec3& position);
Vec3 GetLeadingEdgePosition() const;

void SetChord(Scalar chord);
Scalar GetChord() const;

void SetTwistAngle(Scalar angle_deg);
Scalar GetTwistAngle() const;
Scalar GetTwistAngleRadians() const;
```

## Airfoil

```cpp
void SetAirfoil(std::shared_ptr<Airfoil> airfoil);
std::shared_ptr<const Airfoil> GetAirfoil() const;
bool HasAirfoil() const;
```

## Control Surfaces

```cpp
void AddControlSurface(const ControlSurfaceDefinition& cs);
const std::vector<ControlSurfaceDefinition>& GetControlSurfaces() const;
```

### ControlSurfaceDefinition

```cpp
struct ControlSurfaceDefinition {
    std::string name;
    Scalar hinge_x_fraction;   // hinge location as fraction of chord
    Vec3 hinge_vector;         // hinge axis direction
    Scalar sign;               // deflection sign convention
    Scalar gain;               // deflection gain factor
};
```

## Design Variables

```cpp
void AddDesignVariable(const DesignVariableDefinition& dv);
const std::vector<DesignVariableDefinition>& GetDesignVariables() const;
```

### DesignVariableDefinition

```cpp
struct DesignVariableDefinition {
    std::string name;
    Scalar weight;
};
```
