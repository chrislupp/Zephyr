---
sidebar_position: 5
---

# Control Surfaces

A flap on an externally deformed mesh, with coefficient gradients and a
strip-load adjoint. This is the coupled-simulation shape of the API: a mesh
comes in from somewhere else, a command is applied as a boundary condition, and
derivatives come back for an arbitrary functional of the distributed loads.

```bash
cargo run --example control_surfaces
```

Source: `examples/control_surfaces.rs`.

## The mesh and the command

```rust
use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::Vec3;
use zephyr_vlm::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use zephyr_vlm::solver::{FlightCondition, VortexLatticeSolver};

let mut wing = SurfaceMesh::new("Wing", 4, 6);
wing.y_duplicated = true;
for i in 0..=4 {
    for j in 0..=6 {
        let eta = j as f64 / 6.0;
        let x = i as f64 / 4.0;
        // Representative incoming bending and twist, before actuation.
        wing.set_node(i, j, Vec3::new(x, 3.0 * eta, 0.1 * eta * eta + 0.03 * eta * x))?;
    }
}
wing.controls.push(MeshControlSurface::trailing_edge("flap", 1..6, vec![3.0; 7]));

let mut mesh = MeshInput::default();
mesh.surfaces.push(wing);
mesh.reference.area = 6.0;
mesh.reference.span = 6.0;

let condition = FlightCondition {
    alpha: 4.0_f64.to_radians(),
    control_deflections: vec![("flap".into(), 5.0_f64.to_radians())],
    ..FlightCondition::default()
};

let mut solver = VortexLatticeSolver::new();
solver.set_mesh(&mesh)?;
let solution = solver.solve(&condition)?;
let loads = ForceCalculator::new().compute_with_references(
    &mesh.reference, solver.mesh(), &solution, &condition,
);
println!("CL = {:.6}, CD = {:.6}", loads.coefficients.lift, loads.coefficients.drag);
```

Three things to notice.

The node loop writes a **deformed** shape — a quadratic bend plus a linear
twist. That is the whole point of the direct-mesh interface: camber and twist are
the node positions, so a structural solver's output is the aerodynamic input
without an intermediate parameterization.

`MeshControlSurface::trailing_edge("flap", 1..6, vec![3.0; 7])` puts the hinge
at chordwise node row 3 of 4, across strips 1 through 5, with seven hinge
stations because the mesh has `ny + 1 = 7` spanwise node stations. The hinge
coordinate is **material** — an index into the node grid, not a fraction of the
deformed chord — so it rides along with the deformation.

The deflection goes in the flight condition, not into the mesh. Zephyr rotates
the boundary-condition normals. The mesh handed over is the deformed shape
*before* actuation; if your structural model already rotated those nodes, omit
the command here instead.

`compute_with_references` is the direct-mesh force entry point, taking
`mesh.reference` rather than an `Aircraft`.

## Control gradients

```rust
use zephyr_vlm::gradients::GradientCalculator;

let mut gradients = GradientCalculator::new();
gradients.set_mesh(mesh)?;
for control in gradients.compute_control_gradients(condition)? {
    println!(
        "dCL/d{} = {:.6} per radian",
        control.name, control.by_deflection.coefficients.lift
    );
}
```

Analytical, per radian, for every declared command. Nothing is finite
differenced.

## The strip-load adjoint

This is the part that matters for aeroelasticity. A structural adjoint does not
want $dC_L/d\delta$; it wants the derivative of *its own* functional of the
distributed loads.

```rust
use zephyr_vlm::core::Vec3;
use zephyr_vlm::gradients::{StripLoadAdjoint, StripLoadCotangent};

let count = solver.mesh().strip_count();
let points = vec![Vec3::zeros(); count];

// The functional here is the summed body-z strip force,
// before dynamic pressure.
let mut seed = StripLoadCotangent::zeros(count);
seed.force.fill(Vec3::z());

let adjoint = StripLoadAdjoint::new(solver, condition, &points)?;
println!("{:?}", adjoint.control_cotangents(&[seed.clone()])?);

let mesh_gradients = adjoint.mesh_node_cotangents(mesh, &[seed])?;
println!("input nodes: {}", mesh_gradients[0].surfaces[0].len());
```

A `StripLoadCotangent` holds a force seed and a moment seed per strip. Seeding
every strip's force with $\hat z$ makes the functional the summed vertical force.
A real coupling would seed whatever its own load-transfer operator produces, and
would pass several seeds at once so they share one set of adjoint solves.

Two conventions to get right:

- The seeded quantities are the **raw Kutta–Joukowski sums**, not the normalized
  coefficients `StripLoads` reports. Working in physical units means folding
  dynamic pressure into the seed.
- The moment points — one per strip, `Vec3::zeros()` here — are held **fixed**
  during differentiation. In a real coupling they move with the structure, and
  that motion is a term you differentiate yourself.

`mesh_node_cotangents` returns a cotangent per input node, laid out like
`SurfaceMesh::nodes`. `flight_condition_cotangents` and `control_cotangents`
give the same functional's derivatives with respect to the flight angles and the
commands.

## Restrictions

The adjoint rejects image-plane symmetry and section drag polars rather than
returning an approximate derivative. Mirror with `y_duplicated`, as this example
does. See [Coupled Interfaces](../guide/coupled-interfaces.md).
