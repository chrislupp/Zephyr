---
sidebar_position: 2
---

# Solver Workflow

An `Aircraft` or direct `MeshInput` is converted to a `VortexLattice`.
`VortexLatticeSolver` builds the `3N × N` velocity influence matrix and
`N × N` normalwash matrix, factors the latter, and caches six unit responses.

Changing geometry, symmetry, or Mach invalidates the relevant caches. Changing
only alpha, beta, or rates reuses them. Control deflections rotate affected
panel normals and solve the corresponding control-adjusted linear system.

The returned `SolutionState` owns circulation, six unit responses, and induced
velocity. Pass it with the same mesh and condition to `ForceCalculator`.
Invalid or stale states are rejected by fallible operations.
