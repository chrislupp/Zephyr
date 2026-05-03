---
sidebar_position: 3
---

# Surface

```cpp
#include "pegasus/geometry/surface.h"
```

Represents a lifting surface (wing, tail, fin) composed of spanwise sections.

## Constructors

```cpp
Surface();
Surface(const std::string& name);
```

## Section Management

```cpp
void AddSection(const Section& section);
const std::vector<Section>& GetSections() const;
Index GetNumSections() const;
const Section& GetSection(Index index) const;
```

## Paneling Parameters

```cpp
void SetChordwisePanelCount(int count);
int GetChordwisePanelCount() const;

void SetSpanwisePanelCount(int count);
int GetSpanwisePanelCount() const;

void SetChordwiseSpacing(SpacingType type);
SpacingType GetChordwiseSpacing() const;

void SetSpanwiseSpacing(SpacingType type);
SpacingType GetSpanwiseSpacing() const;
```

## Symmetry and Duplication

```cpp
void SetYDuplicate(Scalar y_plane);
void ClearYDuplicate();
bool IsYDuplicated() const;
```

Calling `SetYDuplicate(0.0)` mirrors the surface about the $y = 0$ plane, creating a symmetric wing from a single half-span definition.

## Coordinate Transformations

```cpp
void SetAngleOffset(Scalar angle_deg);
Scalar GetAngleOffset() const;

void SetTranslation(const Vec3& translation);
Vec3 GetTranslation() const;
```

## Aerodynamic Options

```cpp
void SetHasWake(bool has_wake);
bool HasWake() const;

void SetContributesToLoads(bool contributes);
bool ContributesToLoads() const;
```

## Computed Properties

```cpp
Scalar ComputePlanformArea() const;
Scalar ComputeMeanAerodynamicChord() const;
Scalar ComputeSpan() const;
```
