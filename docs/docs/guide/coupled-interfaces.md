---
sidebar_position: 12
---

# Coupled Interfaces

This sheet is about using Zephyr as the aerodynamic component of something
larger: an aeroelastic solver, a coupled MDO system, or an optimizer that owns
its own geometry parameterization. The library's side of that contract is three
things — a mesh it will panel exactly as given, distributed loads it will
report, and adjoints of arbitrary functionals of those loads.

What Zephyr does **not** do is equally part of the contract. Structural
equilibrium, load transfer, displacement interpolation, total coupled
derivatives, and the fixed-point iteration between disciplines all belong to the
surrounding framework.

## The mesh is the interface

Hand over node coordinates, not design parameters:

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshInput, SurfaceMesh};
use zephyr_vlm::solver::VortexLatticeSolver;

let mut wing = SurfaceMesh::new("Wing", 4, 6);
wing.y_duplicated = true;
for i in 0..=wing.nx {
    for j in 0..=wing.ny {
        let eta = j as f64 / wing.ny as f64;
        let x = i as f64 / wing.nx as f64;
        // Deformed shape from the structural solver, before control actuation.
        wing.set_node(i, j, Vec3::new(x, 3.0 * eta, 0.1 * eta * eta + 0.03 * eta * x))?;
    }
}

let mut mesh = MeshInput::default();
mesh.surfaces.push(wing);
mesh.reference.area = 6.0;
mesh.reference.span = 6.0;

let mut solver = VortexLatticeSolver::new();
solver.set_mesh(&mesh)?;
```

See [Geometry](./geometry.md) for the node layout and for how panels and normals
are built from it. Two points matter here specifically:

- **Camber and twist are the nodes.** There is no airfoil to blend and no twist
  angle to apply; whatever shape the nodes describe is the shape that is solved.
  That is what makes a deformed mesh mean what the structural solver intended.
- **Pass the deformed mesh before control actuation**, and put the deflection in
  `FlightCondition::control_deflections`. Zephyr rotates the boundary-condition
  normals; it does not move nodes or the wake. If your structural model already
  physically rotated the control-surface nodes, omit the command here — see
  [Control Surfaces](./control-surfaces.md).

## Distributed loads out

`AerodynamicResults::surface_loads` carries a `StripLoads` for every strip of
every generated surface, mirrors included, in lattice order. That is the load
distribution to transfer. [Loads and Coefficients](./loads.md) lists the fields
and the normalization caveat: strip coefficients divide by the **aircraft**
reference area and chord, not by the strip's own.

## `StripLoadAdjoint`: derivatives of arbitrary load functionals

Coefficient gradients are not enough for a coupled system. What a structural
adjoint needs is the derivative of *its own* functional of the distributed
loads — a weighted stress, a displacement, an aggregated margin — with respect
to the mesh, the flight condition, and the control commands.

`StripLoadAdjoint` takes the seed of that functional and pulls it back:

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::gradients::{StripLoadAdjoint, StripLoadCotangent};

let count = solver.mesh().strip_count();
let moment_points = vec![Vec3::zeros(); count];

// One functional: the summed body-z strip force.
let mut seed = StripLoadCotangent::zeros(count);
seed.force.fill(Vec3::z());

let adjoint = StripLoadAdjoint::new(solver, condition, &moment_points)?;

let by_controls = adjoint.control_cotangents(&[seed.clone()])?;
let by_flight = adjoint.flight_condition_cotangents(&[seed.clone()])?;
let by_nodes = adjoint.mesh_node_cotangents(mesh, &[seed])?;

println!("surfaces: {}", by_nodes[0].surfaces.len());
println!("nodes on surface 0: {}", by_nodes[0].surfaces[0].len());
```

### Seeds

A `StripLoadCotangent` is one functional's seed: a `Vec3` per strip for the
resultant force and a `Vec3` per strip for the resultant moment.
`StripLoadCotangent::zeros(strip_count)` starts from nothing. Every method takes
a **slice** of them, so several functionals share one set of adjoint solves.

The convention to get right: the seeded quantities are the **raw
Kutta–Joukowski sums** — the sum of a strip's panel forces, and their moment
about that strip's moment point — not the area- and chord-normalized
coefficients `StripLoads` reports. A caller working in physical units folds its
dynamic pressure into the cotangents.

### Moment points

`StripLoadAdjoint::new` takes one moment point per strip, which is where each
strip's moment is taken. There must be exactly `strip_count` of them.

Those points are held **fixed** during differentiation. If your functional's
moment points move with the structure — and in an aeroelastic coupling they
usually do — that motion is a separate term you differentiate yourself.

### What comes back

| Method | Returns |
| --- | --- |
| `mesh_node_cotangents(&mesh, &seeds)` | one `MeshNodeCotangent` per seed, each with a `Vec3` per node per input surface, laid out like `SurfaceMesh::nodes` |
| `flight_condition_cotangents(&seeds)` | one `TotalDerivative` per seed, with `by_alpha` and `by_beta` |
| `control_cotangents(&seeds)` | per seed, a list of `(name, derivative)` pairs per radian, sorted like the lattice's control list |

The mesh-node pullback requires the `MeshInput` to be the one that produced the
solved lattice; a mismatch is `ZephyrError::InvalidMesh`, not a wrong answer.
A cotangent whose length does not match the strip count is
`ZephyrError::InvalidGeometry`. A stale solution is `ZephyrError::InvalidSolution`.

### Restrictions

The same two that apply to geometry adjoints, for the same reasons, and both are
errors rather than silent approximations:

- **Explicit mirroring only.** A lattice using `SymmetryConfig` image planes is
  rejected. Duplicate the wing with `y_duplicated` instead.
- **No section drag polars.** A viscous strip's drag depends on its own lift
  coefficient and therefore on the strip area; those tangents do not exist yet.

Control cotangents include the induced velocity and use the solved, deflected
AIC, so they are valid at nonzero deflection.

## Trim inside the coupling

An external framework can either supply trim commands itself or let Zephyr solve
a local trim and hand back derivatives of the converged map. The second is
`MeshTrimSolver` — see [Trim](./trim.md). The key properties for a coupled
solver:

- External loads are **fixed inside** each Zephyr solve. A load that depends on
  a solved operating variable has to be connected through your model and
  converged there.
- All function gradients include the trim response, including derivatives with
  respect to deformed node coordinates, targets, external loads, and dynamic
  pressure.
- Trim needs a nonsingular local Jacobian. Ineffective controls or infeasible
  bounds are reported as failures, not approximated.

## Over a wire

For frameworks that talk gRPC, the `philote` feature exposes the same thing as
Philote-MDO disciplines — a mesh discipline and a VLM discipline, both explicit,
with analytical partials and optional internal trim. See
[MDO Integration](./mdo-integration.md).

Distributed loads are **not** currently transported through the Philote
discipline; its outputs are integrated coefficients and resultants. Strip-level
coupling goes through the native API above.

## Example

`cargo run --example control_surfaces` builds a deformed mesh with a flap,
reports coefficients and control gradients, and then takes a strip-load adjoint
against a summed-force functional. See
[Control Surfaces (example)](../examples/control-surfaces.md).
