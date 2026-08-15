---
sidebar_position: 9
---

# EigenmodeAnalysis

`EigenmodeAnalysis` accepts aircraft references, `MassProperties`,
`StabilityDerivatives`, and `FlightState`. It builds longitudinal/lateral
4-state or full 12-state matrices and returns normalized complex eigenpairs
classified as short period, phugoid, Dutch roll, roll subsidence, or spiral.
