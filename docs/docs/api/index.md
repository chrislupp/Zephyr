---
sidebar_position: 1
---

# API Reference

VLM provides a single facade header `vlm/vlm.h` that includes all public API components. For finer-grained control, include individual headers directly.

All types are in the `vlm` namespace.

## Core Types

Defined in `vlm/core/types.h`:

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
| [`Aircraft`](aircraft.md) | `vlm/geometry/aircraft.h` | Top-level aircraft configuration |
| [`Surface`](surface.md) | `vlm/geometry/surface.h` | Lifting surface definition |
| [`Section`](section.md) | `vlm/geometry/section.h` | Wing section at a spanwise station |
| [`VortexLatticeSolver`](vortex-lattice-solver.md) | `vlm/solver/vortex_lattice_solver.h` | Core VLM solver |
| [`ForceCalculator`](force-calculator.md) | `vlm/aero/force_calculator.h` | Aerodynamic force computation |
| [`StabilityDerivativeCalculator`](stability-derivatives.md) | `vlm/aero/stability_derivatives.h` | Stability derivative computation |
| [`TrimSolver`](trim-solver.md) | `vlm/trim/trim_solver.h` | Newton-Raphson trim solver |
| [`EigenmodeAnalysis`](eigenmode-analysis.md) | `vlm/dynamics/eigenmode_analysis.h` | Dynamic stability eigenmode analysis |
| [`AvlReader`](avl-reader.md) | `vlm/io/avl_reader.h` | AVL geometry file parser |
| [`GradientCalculator`](gradient-calculator.md) | `vlm/ad/gradient_calculator.h` | Analytical gradient computation |
