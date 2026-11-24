---
sidebar_position: 6
---

# Mesh Trim

A deformed wing and tail trimmed to vertical-force and pitching-moment
equilibrium by solving for angle of attack and elevator, followed by exact
derivatives of the converged mapping.

```bash
cargo run --example mesh_trim
```

Source: `examples/mesh_trim.rs`.

## The geometry and the problem

```rust
use zephyr_vlm::{
    core::Vec3,
    geometry::{MeshControlSurface, MeshInput, SurfaceMesh},
    solver::FlightCondition,
    trim::{AeroQuantity, MeshTrimProblem, MeshTrimSolver, TrimBalance, TrimVariable},
};

let mut mesh = MeshInput::default();
for (name, xle, span, chord) in [("Wing", 0.0, 3.0, 1.0), ("Tail", 3.0, 1.3, 0.65)] {
    let mut surface = SurfaceMesh::new(name, 2, 3);
    surface.y_duplicated = true;
    for i in 0..=2 {
        for j in 0..=3 {
            let eta = j as f64 / 3.0;
            surface.set_node(
                i, j,
                Vec3::new(xle + chord * i as f64 / 2.0, span * eta, 0.08 * eta * eta),
            )?;
        }
    }
    if name == "Tail" {
        surface.controls.push(MeshControlSurface::trailing_edge(
            "elevator", 0..3, vec![1.2; 4],
        ));
    }
    mesh.surfaces.push(surface);
}
mesh.reference.area = 7.0;
mesh.reference.chord = 0.9;
mesh.reference.span = 6.0;
mesh.reference.moment_reference_point = Vec3::new(0.3, 0.0, 0.0);

let mut alpha = TrimVariable::new("alpha", 0.03);
alpha.bounds = [-0.3, 0.3];
let mut elevator = TrimVariable::new("control_elevator", 0.0);
elevator.bounds = [-0.4, 0.4];

let problem = MeshTrimProblem {
    variables: vec![alpha, elevator],
    balances: vec![
        TrimBalance::new(AeroQuantity::Fz, 1000.0),
        TrimBalance::new(AeroQuantity::My, 1000.0),
    ],
    ..MeshTrimProblem::default()
};

let condition = FlightCondition {
    mach: 0.25,
    dynamic_pressure: 650.0,
    ..FlightCondition::default()
};
// Prescribed external body-z load. This case enforces Fz and My only.
let external = [0.0, 0.0, -1200.0, 0.0, 0.0, 0.0];

let solution = MeshTrimSolver.solve(&mesh, &condition, &problem, &[0.0, 0.0], &external)?;
```

Both surfaces are deformed half models with a quadratic dihedral term; the tail
carries the elevator. Two unknowns, two equations — the system must be square,
with exactly one independent variable per unique balance quantity.

The balances are **dimensional**: each row is aerodynamic load plus prescribed
external load minus target. With an external $F_z$ of $-1200$ N and a zero
target, the trim drives the aerodynamic $F_z$ to $+1200$ N. The `1000.0`
arguments are residual scales in the quantity's own units, not targets.

External loads use body axes and mathematical right-hand moments about the
configured reference point, so `Mx` and `Mz` carry none of the conventional
`Cl`/`Cn` sign flips. Getting that wrong is the most likely way to trim to
something you did not intend.

Speed, Mach, dynamic pressure, and the rates stay prescribed. This case enforces
two of the six equilibrium equations, not aircraft equilibrium.

## Reading the result

```rust
println!("alpha = {:.6} rad, elevator = {:.6} rad",
    solution.states[0], solution.states[1]);
println!("CL = {:.6}, CD = {:.6}", solution.values[0], solution.values[1]);
println!("physical residuals = {:?}", solution.diagnostics.residuals);
```

`values` is all twelve `AeroQuantity` values in order — coefficients first, then
dimensional loads. `diagnostics.residuals` are **unscaled**, in physical units,
which is what you want when deciding whether a converged solve is actually tight
enough. `diagnostics.singular_value_ratio` is the conditioning of the scaled
Jacobian; a very small value means the controls are nearly redundant for these
balances.

A failed trim is an error, not a status. It never returns an approximate
derivative.

## Derivatives through trim

```rust
let partials = solution.compute_partials(
    mesh, &[AeroQuantity::CL, AeroQuantity::CD], true,
)?;

println!("dCD/d(dynamic_pressure) = {:.9} /Pa",
    partials[&("CD".into(), "dynamic_pressure".into())][0]);
println!("dalpha/d(external_Fz) = {:.9} rad/N",
    partials[&("alpha".into(), "external_Fz".into())][0]);
```

The first of those two is the one worth staring at. $C_D$ has **no** direct
dynamic-pressure derivative — coefficients are dimensionless. It is nonzero here
because raising the dynamic pressure changes the trimmed angle of attack through
the dimensional $F_z$ balance, and the drag follows. Every block
`compute_partials` returns includes that trim response.

Solved trim states are outputs too, which is why `dalpha/d(external_Fz)` is
available alongside the functions you asked for.

The third argument controls whether mesh-node blocks are computed. Pass `false`
inside an iteration that holds geometry fixed and you skip the mesh adjoint
work entirely.

`compute_partials_adjoint` computes the same thing with transpose solves — the
cheaper direction when there are many inputs and few functions.

## Next

- [Trim](../guide/trim.md) for control mixers, the numerics, and what makes the
  mapping differentiable.
- [MDO Integration](../guide/mdo-integration.md) to run this same trim inside a
  gRPC discipline.
- [Philote Optimization](./philote-optimization.md) for an optimizer driving it.
