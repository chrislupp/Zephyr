---
sidebar_position: 3
---

# Stability Analysis

`StabilityDerivativeCalculator` combines cached unit responses into
derivatives with respect to alpha, beta, normalized roll/pitch/yaw rates, and
named controls. Results use grouped `Coefficients`, avoiding the ambiguous
`CL`/`Cl` naming found in the former C++ API.

`TrimSolver` assembles these derivatives into a Newton Jacobian. Nonconvergence
is returned as `TrimResult` with status, iteration count, residual, and a
diagnostic.

`EigenmodeAnalysis` constructs longitudinal and lateral 4-state systems or a
full 12-state rigid-body system. It returns normalized complex eigenvectors,
deterministically sorted and classified modes, and the exact system matrix.
