---
sidebar_position: 14
---

# StripLoadAdjoint

`gradients::StripLoadAdjoint` differentiates arbitrary functionals of the
**distributed** strip forces and moments — the interface a coupled aeroelastic
adjoint needs, where coefficient gradients are not enough.

The guide sheet is [Coupled Interfaces](../guide/coupled-interfaces.md); the
worked example is [Control Surfaces](../examples/control-surfaces.md).

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::gradients::{StripLoadAdjoint, StripLoadCotangent};

let count = solver.mesh().strip_count();
let moment_points = vec![Vec3::zeros(); count];

let mut seed = StripLoadCotangent::zeros(count);
seed.force.fill(Vec3::z());

let adjoint = StripLoadAdjoint::new(solver, condition, &moment_points)?;
let by_nodes = adjoint.mesh_node_cotangents(mesh, &[seed])?;
```

## Construction

```text
StripLoadAdjoint::new(&solver, &condition, &moment_points) -> Result<Self>
```

`moment_points` is one point per strip, where that strip's moment is taken.
There must be exactly `strip_count` of them.

Those points are held **fixed** during differentiation. In an aeroelastic
coupling they usually move with the structure, and that motion is a separate
term you differentiate yourself.

Errors: `ZephyrError::InvalidSolution` for a stale solution,
`ZephyrError::InvalidGeometry` for a moment-point count that does not match the
lattice or for an unsupported configuration below.

`mesh()` returns the lattice the adjoint is built on.

## `StripLoadCotangent`

One functional's seed: `force` and `moment`, each a `Vec3` per strip.
`zeros(strip_count)` starts from nothing.

The seeded quantities are the **raw Kutta–Joukowski sums** — the sum of a
strip's panel forces, and their moment about that strip's moment point — not the
area- and chord-normalized coefficients `StripLoads` reports. A caller working
in physical units folds its dynamic pressure into the cotangents.

Every method takes a **slice** of seeds, so several functionals share one set of
adjoint solves. Pass them together.

## Pullbacks

| Method | Returns |
| --- | --- |
| `mesh_node_cotangents(&MeshInput, &seeds)` | one `MeshNodeCotangent` per seed: a `Vec3` per node per input surface, laid out like `SurfaceMesh::nodes` |
| `flight_condition_cotangents(&seeds)` | one `TotalDerivative` per seed, with `by_alpha` and `by_beta` |
| `control_cotangents(&seeds)` | per seed, `(name, derivative)` pairs per radian, sorted like the lattice's control list |

`mesh_node_cotangents` requires the `MeshInput` to be the one that produced the
solved lattice; a mismatch is `ZephyrError::InvalidMesh`, not a wrong answer. A
cotangent whose length does not match the strip count is
`ZephyrError::InvalidGeometry`.

Control cotangents include the induced velocity and use the solved, deflected
AIC, so they are valid at nonzero deflection.

## Restrictions

The same two as the geometry adjoints, and both are errors rather than silent
approximations:

- **Explicit mirroring only.** A lattice using `SymmetryConfig` image planes is
  rejected; duplicate the wing with `y_duplicated`.
- **No section drag polars.** A viscous strip's drag depends on its own lift
  coefficient and therefore on the strip area; those tangents do not exist yet.

[`StripLoadAdjoint` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/gradients/struct.StripLoadAdjoint.html)
