---
sidebar_position: 6
---

# VortexLatticeSolver

`solver::VortexLatticeSolver` holds a lattice, two influence operators, a
factorization, six unit responses, and the last solution. The cache behavior is
the thing to understand; see
[Solver Workflow](../guide/solver-workflow.md) for the full invalidation table.

```rust
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(aircraft)?;
let solution = solver.solve(condition)?;
```

## Installing geometry

| Method | Effect |
| --- | --- |
| `set_aircraft(&Aircraft)` | mesh a section model; adopts its references, moment point, and Mach |
| `set_mesh(&MeshInput)` | install a direct mesh; adopts its references and moment point |
| `set_reference_point(Vec3)` | move the moment and rate-rotation point; keeps the factorization, rebuilds the unit responses |
| `mutate_mesh(closure)` | scoped in-place lattice mutation, invalidating everything dependent |

Each of the first two discards everything cached from the other. `mutate_mesh`
is the only supported way to change a lattice in place; mutating through any
other route leaves caches claiming a validity they do not have.

## Solving

| Method | Effect |
| --- | --- |
| `solve(&FlightCondition)` | the whole sequence, doing only what is stale |
| `rebuild_aic()` | assemble both operators and factor, explicitly |
| `solve_unit_responses()` | solve the six unit responses, explicitly |

`solve` rejects non-finite angles, rates, Mach or dynamic pressure, negative
Mach or dynamic pressure, unknown control names, and any Mach at or above 1 —
Prandtl–Glauert has no meaning there.

The two explicit methods exist so you can pay the assembly and factorization
cost at a predictable moment rather than inside the first `solve`.

## Symmetry and accelerations

| Method | Effect |
| --- | --- |
| `set_symmetry(SymmetryConfig)` | image planes: influence without unknowns. **Unsupported by the geometry and strip-load adjoints.** |
| `set_mirror_reuse(bool)` | derive mirrored influence entries by reflection. Bit-identical results. Off by default. |
| `set_symmetric_factorization(bool)` | factor two half-size systems when the matrix is exactly mirror-symmetric. Agrees to rounding, not bit for bit. Off by default. |
| `symmetric_factorization_applies()` | whether the last solve actually used the half-size form |

See [Performance](../guide/performance.md).

## Controls

`control_circulation_derivative(condition, name)` returns the analytical
circulation derivative for one named command at the current condition — one
triangular solve, not a re-solve. It is what the control derivatives are built
from.

A deflected solve builds and factors a second matrix from the rotated normals
and caches it **against those normals**, so sweeping alpha at a fixed deflection
reuses it and returning to a previous deflection reuses it again.

## Inspection

| Accessor | Gives |
| --- | --- |
| `mesh()` | the `VortexLattice`: panels, strips, surface ranges, control names |
| `solution()` | the last `SolutionState` |
| `aic()` | the undeflected normalwash matrix |
| `solution_aic()`, `linearized_aic()` | the matrices the last solve and linearization used |
| `velocity_matrix()` | the `3N × N` force-velocity operator, self-bound leg excluded |
| `control_velocity_matrix()` | the velocity operator at control points |
| `reference_quantities()`, `reference_point()` | the normalization in force |
| `current_mach()` | the Mach the matrices were assembled at |
| `is_aic_valid()`, `are_unit_responses_valid()` | cache state |

`freestream_direction(&condition)` is an associated function returning the unit
freestream vector, for code that needs to agree with the solver's convention.

## Two operators, not one

The normalwash matrix includes **every** bound vortex at the control points; the
velocity operator excludes each panel's **own** bound leg and is evaluated at
the bound-force points. Including the self-bound leg is correct for the boundary
condition and wrong for the Kutta–Joukowski force. Both are assembled in the
same pass, and the adjoints differentiate both.

## `SolutionState`

The combined `circulation`, all six unit responses
(`circulation_per_u` … `circulation_per_r`), the `N × 3` `induced_velocity` in
body axes, and `is_valid`. Anything that would silently consume a stale state
rejects it with `ZephyrError::InvalidSolution` instead.

[`VortexLatticeSolver` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/solver/struct.VortexLatticeSolver.html)
