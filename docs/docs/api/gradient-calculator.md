---
sidebar_position: 11
---

# GradientCalculator

Set an `Aircraft` or `MeshInput`, choose outputs with `OutputSelection`, then
request flight-condition, section-geometry, or mesh-node gradients. Scalar
objectives use `AeroOutput::{Lift, Drag, Side, Roll, Pitch, Yaw,
InducedDrag}`.

Low-level APIs expose safe Biot–Savart VJPs, AIC row/column derivatives,
near-field force derivatives, Trefftz derivatives, Prandtl–Glauert/Mach
derivatives, and batched transpose solves backed by the solver's cached LU
factorization.
