---
sidebar_position: 4
---

# Solver Workflow

`VortexLatticeSolver` is a stateful object that holds a lattice, two influence
operators, a factorization, six unit responses, and the last solution. Nearly
all of the performance available in this library comes from understanding what
of that survives a change and what does not.

The short version: **keep one solver alive for as long as the geometry and Mach
hold still.**

## The sequence

```rust
let mut solver = VortexLatticeSolver::new();
solver.set_aircraft(aircraft)?;        // generate the lattice
let solution = solver.solve(condition)?; // assemble, factor, solve
let results = ForceCalculator::new()
    .compute_for_aircraft(aircraft, solver.mesh(), &solution, condition);
```

`set_aircraft` meshes a section model; `set_mesh` takes a `MeshInput` directly.
Each discards everything cached from the other, and each adopts the reference
quantities and moment reference point of what it was given — `set_aircraft` also
adopts the aircraft's Mach number as the solver's current Mach.

`solve` does as little as it can get away with:

1. Validate the condition. Non-finite angles, rates, Mach or dynamic pressure
   are rejected, as are negative Mach and dynamic pressure, unknown control
   names, and any Mach at or above 1 — Prandtl–Glauert has no meaning there, and
   returning the incompressible answer silently would hand back plausible
   numbers for a condition the method cannot model.
2. If Mach changed, invalidate the influence matrices.
3. If the influence matrices are stale, assemble both of them and factor the
   normalwash matrix.
4. If the unit responses are stale, solve the six of them in one batched
   back-substitution.
5. Combine them for this condition, or — if controls are deflected — solve the
   control-adjusted system.

`rebuild_aic` and `solve_unit_responses` are public if you want to pay those
costs at a predictable moment rather than inside the first `solve`.

## The six unit responses

The right-hand side is linear in the onset velocity, so the circulation is a
linear combination of six precomputed solutions: three for unit freestream
components $u$, $v$, $w$, and three for unit body rates $p$, $q$, $r$.

$$
\boldsymbol{\Gamma} = \hat V_x \boldsymbol{\Gamma}_u + \hat V_y \boldsymbol{\Gamma}_v
+ \hat V_z \boldsymbol{\Gamma}_w + p\,\boldsymbol{\Gamma}_p + q\,\boldsymbol{\Gamma}_q
+ r\,\boldsymbol{\Gamma}_r
$$

Once they are in hand, a change of alpha, beta, or any rate costs six scaled
vector adds — no triangular solve at all. They are also exactly what the
stability derivatives are assembled from, which is why
`StabilityDerivativeCalculator::compute` needs no extra solves.

`SolutionState` carries the combined circulation, all six responses, the induced
velocity, and an `is_valid` flag. Anything that would silently consume a stale
state rejects it instead.

## What invalidates what

| Change | Lattice | Influence matrices and LU | Unit responses | Deflected LU |
| --- | --- | --- | --- | --- |
| `alpha`, `beta`, body rates | kept | kept | kept | kept |
| `dynamic_pressure` | kept | kept | kept | kept |
| `mach` | kept | **rebuilt** | **rebuilt** | **rebuilt** |
| control deflections | kept | kept | kept | **refactored** |
| `set_symmetry` | kept | **rebuilt** | **rebuilt** | **rebuilt** |
| `set_reference_point` | kept | kept | **rebuilt** | kept |
| `mutate_mesh` | mutated | **rebuilt** | **rebuilt** | **rebuilt** |
| `set_aircraft` / `set_mesh` | **regenerated** | **rebuilt** | **rebuilt** | **rebuilt** |

Dynamic pressure never enters the solve: circulation is computed against a unit
freestream, and dynamic pressure appears only where a dimensional load is
formed. That is why an alpha sweep at several dynamic pressures costs exactly
one factorization.

`mutate_mesh` takes a closure over the `VortexLattice` and invalidates
everything dependent on it afterwards. It is the supported way to nudge a
lattice in place; mutating through any other route leaves caches claiming
validity they do not have.

## Deflected controls

With a nonempty deflection list, the panel normals are rotated and a second
matrix is built and factored from the rotated normals. That factorization is
cached against the rotated normals themselves, so sweeping alpha at a fixed
deflection reuses it, and returning to the same deflection later reuses it
again. Changing a deflection refactors.

A deflected solve therefore costs about one extra factorization per distinct
deflection state — which is why a trim iteration over alpha and elevator is
dominated by the elevator steps, not the alpha steps.

## Symmetry and the two accelerations

`set_symmetry` installs image planes (`SymmetryConfig`) that contribute
influence without adding unknowns. It is the cheapest form of symmetry and the
least capable: the geometry and strip-load adjoints do not support it.

For lattices built from explicitly duplicated surfaces, two opt-in
accelerations are available, both off by default:

- `set_mirror_reuse(true)` derives the mirrored half of both influence operators
  by reflection instead of evaluating the kernel for it. Results are
  **bit-identical** to having it off, because the kernel is reflection
  equivariant in floating point by construction.
- `set_symmetric_factorization(true)` factors two half-size systems instead of
  one full one when the system matrix is exactly mirror-symmetric. Two
  half-size factorizations cost a quarter of the full one. Results agree with
  the full factorization to the rounding of the elimination, not bit for bit.

The symmetric form does not hold for antisymmetric deflections such as an
aileron, for partially duplicated geometry, or for a half-wing model; those
cases fall back to the full LU silently.
`symmetric_factorization_applies()` reports which happened for the last solve.

[Performance](./performance.md) covers both in more detail, along with the
`parallel` feature and the half-wing mode.

## Inspecting the solver

| Accessor | Gives |
| --- | --- |
| `mesh()` | the `VortexLattice`: panels, strips, surface ranges, control names |
| `aic()` | the undeflected normalwash matrix |
| `solution_aic()` / `linearized_aic()` | the matrix the last solve and the last linearization actually used |
| `velocity_matrix()` | the `3N × N` force-velocity operator |
| `control_velocity_matrix()` | the velocity operator at control points |
| `solution()` | the last `SolutionState` |
| `is_aic_valid()` / `are_unit_responses_valid()` | cache state |
| `reference_quantities()` / `reference_point()` | normalization in force |
| `current_mach()` | the Mach the matrices were assembled at |

These are the entry points for building something the library does not provide.
The adjoint machinery uses nothing else.

## When not to use this object

`VortexLatticeSolver` always computes a complete solution: both influence
operators, all six unit responses, the full circulation. If you only need two
coefficients and one derivative, that is more work than you asked for. The
request-driven interface in [Selective Analysis](./selective-analysis.md) skips
the stages a request does not reach, and can solve only one semispan of a
symmetric model.
