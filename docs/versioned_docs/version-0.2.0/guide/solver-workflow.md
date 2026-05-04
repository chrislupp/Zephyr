---
sidebar_position: 2
---

# Solver Workflow

This guide explains the complete VLM analysis pipeline from geometry definition to results extraction.

## Geometry Hierarchy

VLM uses a three-level geometry hierarchy:

```
Aircraft
├── Surface (wing, tail, fin, ...)
│   ├── Section (root, break, tip)
│   └── Section
├── Surface
│   ├── Section
│   └── Section
└── ...
```

### Aircraft

The `Aircraft` class is the top-level container. It holds all lifting surfaces, reference quantities for non-dimensionalization, and global settings.

```cpp
Aircraft aircraft("My Aircraft");
aircraft.SetReferenceArea(30.0);    // m^2
aircraft.SetReferenceChord(2.0);    // m
aircraft.SetReferenceSpan(15.0);    // m
aircraft.SetMomentReferencePoint(Vec3(0.5, 0.0, 0.0));
```

### Surface

Each `Surface` represents a lifting surface with specified paneling parameters:

```cpp
Surface wing("Main Wing");
wing.SetChordwisePanelCount(8);
wing.SetSpanwisePanelCount(30);
wing.SetChordwiseSpacing(SpacingType::kCosine);
wing.SetSpanwiseSpacing(SpacingType::kCosine);
wing.SetYDuplicate(0.0);  // symmetric about XZ plane
```

### Section

`Section` objects define the spanwise stations of each surface:

```cpp
Section root(Vec3(0.0, 0.0, 0.0), 3.0, 2.0);  // LE position, chord, twist (deg)
Section tip(Vec3(1.0, 7.5, 0.5), 1.5, -1.0);
wing.AddSection(root);
wing.AddSection(tip);
```

## Solver Setup

The `VortexLatticeSolver` builds the mesh and aerodynamic influence coefficient (AIC) matrix:

```cpp
VortexLatticeSolver solver;
solver.SetAircraft(aircraft);
```

The AIC matrix is built lazily on the first solve. To explicitly trigger it:

```cpp
solver.RebuildAIC();
```

## Flight Conditions

The `FlightCondition` struct specifies the flow state:

```cpp
FlightCondition condition;
condition.alpha = 5.0 * kDegreesToRadians;  // angle of attack
condition.beta = 0.0;                        // sideslip angle
condition.roll_rate = 0.0;                   // non-dimensional roll rate (pb/2V)
condition.pitch_rate = 0.0;                  // non-dimensional pitch rate (qc/2V)
condition.yaw_rate = 0.0;                    // non-dimensional yaw rate (rb/2V)
```

## Solving

```cpp
SolutionState solution = solver.Solve(condition);
```

The solution contains the circulation distribution across all panels.

## Force Computation

```cpp
ForceCalculator calculator;
AerodynamicResults results = calculator.ComputeForces(
    aircraft, solver.GetMesh(), solution, condition);
```

The `AerodynamicResults` struct contains:
- Lift, drag, side force coefficients: `CL`, `CD`, `CY`
- Moment coefficients: `Cl` (roll), `Cm` (pitch), `Cn` (yaw)
- Induced drag from Trefftz plane: `CDi`
- Per-surface force breakdowns

## Control Surfaces

Control surfaces are defined on sections:

```cpp
ControlSurfaceDefinition elevator;
elevator.name = "elevator";
elevator.hinge_x_fraction = 0.7;  // hinge at 70% chord
elevator.hinge_vector = Vec3(0.0, 1.0, 0.0);
elevator.sign = 1.0;
elevator.gain = 1.0;

section.AddControlSurface(elevator);
```

To deflect control surfaces in a flight condition:

```cpp
condition.control_deflections["elevator"] = 5.0 * kDegreesToRadians;
```
