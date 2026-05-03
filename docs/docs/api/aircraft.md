---
sidebar_position: 2
---

# Aircraft

```cpp
#include "vlm/geometry/aircraft.h"
```

Top-level geometry container representing a complete aircraft configuration.

## Constructors

```cpp
Aircraft();
Aircraft(const std::string& title);
```

## Reference Quantities

```cpp
void SetReferenceArea(Scalar area);
Scalar GetReferenceArea() const;

void SetReferenceChord(Scalar chord);
Scalar GetReferenceChord() const;

void SetReferenceSpan(Scalar span);
Scalar GetReferenceSpan() const;

void SetMomentReferencePoint(const Vec3& point);
Vec3 GetMomentReferencePoint() const;
```

## Flow Settings

```cpp
void SetMachNumber(Scalar mach);
Scalar GetMachNumber() const;

void SetProfileDragCoefficient(Scalar cd0);
Scalar GetProfileDragCoefficient() const;
```

## Surface Management

```cpp
void AddSurface(const Surface& surface);
const std::vector<Surface>& GetSurfaces() const;
Index GetNumSurfaces() const;
const Surface& GetSurface(Index index) const;
const Surface* FindSurfaceByName(const std::string& name) const;
```

## Computed Properties

```cpp
Index ComputeTotalVortexCount() const;
Index ComputeTotalStripCount() const;
Scalar ComputeTotalPlanformArea() const;
Scalar ComputeTotalSpan() const;
```

## Control Surfaces

```cpp
std::vector<std::string> GetControlSurfaceNames() const;
Index GetNumControlSurfaces() const;
```
