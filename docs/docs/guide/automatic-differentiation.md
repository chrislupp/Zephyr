---
sidebar_position: 4
---

# Analytical and Adjoint Gradients

Zephyr provides safe analytical VJPs for the regularized Biot–Savart kernel,
AIC rows and columns, near-field forces, and Trefftz induced drag.
`TotalDerivative` solves the transposed aerodynamic system for an `AeroOutput`.

The alpha derivative contains both terms:

```text
dJ/dalpha = explicit fixed-circulation force/axis term
          + lambdaᵀ db/dalpha
```

The explicit term includes freestream dependence and the
body-to-stability-axis rotation. This completes the derivative that was
disabled in Pegasus 0.2.0.

`GradientCalculator` supports flight-condition, aircraft-section, and direct
mesh-node inputs. Angle, rate, and Mach derivatives use analytical influence
and force kernels. Geometry evaluation performs one primal solve, reuses the
cached transpose factorization for a batched adjoint, and contracts each
section or node direction without resolving the aerodynamic system.

`OutputSelection` leaves disabled results exactly zero and avoids their
adjoint and quadratic velocity-kernel work. The deprecated
finite-difference-step setter is a compatibility no-op.
