---
sidebar_position: 1
---

# API Reference

VLM provides a single facade header `pegasus/vlm.h` that includes all public API components. For finer-grained control, include individual headers directly.

All types are in the `vlm` namespace.

## Core Types

Defined in `pegasus/core/types.h`:

| Type | Definition |
|------|-----------|
| `Scalar` | `double` |
| `Vec2`, `Vec3`, `Vec4` | Fixed-size Eigen vectors |
| `Mat2`, `Mat3`, `Mat4` | Fixed-size Eigen matrices |
| `VecX`, `MatX` | Dynamic-size Eigen types |
| `Index` | `Eigen::Index` |

## Enums

| Enum | Values | Description |
|------|--------|-------------|
| `SymmetryType` | `kNone`, `kSymmetric`, `kAntisymmetric` | Flow symmetry type |
| `SpacingType` | `kNone`, `kUniform`, `kCosine`, `kSine`, `kBlend` | Panel spacing distribution |

## Constants

| Constant | Description |
|----------|-------------|
| `kPi` | $\pi$ |
| `kDegreesToRadians` | $\pi / 180$ |
| `kRadiansToDegrees` | $180 / \pi$ |
| `kEpsilon` | Machine epsilon |
| `kDefaultTolerance` | Default solver tolerance |

## Class Overview

| Class | Header | Description |
|-------|--------|-------------|
| [`Aircraft`](aircraft.md) | `pegasus/geometry/aircraft.h` | Top-level aircraft configuration |
| [`Surface`](surface.md) | `pegasus/geometry/surface.h` | Lifting surface definition |
| [`Section`](section.md) | `pegasus/geometry/section.h` | Wing section at a spanwise station |
| [`VortexLatticeSolver`](vortex-lattice-solver.md) | `pegasus/solver/vortex_lattice_solver.h` | Core VLM solver |
| [`ForceCalculator`](force-calculator.md) | `pegasus/aero/force_calculator.h` | Aerodynamic force computation |
| [`StabilityDerivativeCalculator`](stability-derivatives.md) | `pegasus/aero/stability_derivatives.h` | Stability derivative computation |
| [`TrimSolver`](trim-solver.md) | `pegasus/trim/trim_solver.h` | Newton-Raphson trim solver |
| [`EigenmodeAnalysis`](eigenmode-analysis.md) | `pegasus/dynamics/eigenmode_analysis.h` | Dynamic stability eigenmode analysis |
| [`AvlReader`](avl-reader.md) | `pegasus/io/avl_reader.h` | AVL geometry file parser |
| [`GradientCalculator`](gradient-calculator.md) | `pegasus/gradients/gradient_calculator.h` | Analytical gradient computation |
