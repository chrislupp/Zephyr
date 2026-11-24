---
sidebar_position: 10
---

# Analytical and Adjoint Gradients

Zephyr differentiates itself by hand. There is no AD framework, no dual-number
pass, and — apart from the tests, where finite differences are the reference the
analytical kernels are checked against — no finite differencing anywhere in the
library. Every derivative is a closed-form vector–Jacobian product through the
same kernels the primal solve uses.

Two consequences follow, and they are the reason to prefer this interface:

- **Cost.** A scalar output's gradient with respect to any number of inputs
  costs one transposed solve, reusing the factorization the primal already
  produced. A hundred mesh-node coordinates cost one adjoint, not a hundred
  re-solves.
- **Accuracy.** There is no step size to pick, and no cancellation at small
  steps.

## The decomposition

The circulation solve is the only implicit step. Everything else is explicit.
So for a scalar output $J$ and an input $x$,

$$
\frac{dJ}{dx} = \underbrace{\frac{\partial J}{\partial x}\bigg|_{\boldsymbol{\Gamma}}}_{\text{explicit}}
+ \underbrace{\boldsymbol{\lambda}^{\!\top}\!\left(\frac{\partial \mathbf b}{\partial x} - \frac{\partial A}{\partial x}\boldsymbol{\Gamma}\right)}_{\text{implicit}},
\qquad
A^{\!\top}\boldsymbol{\lambda} = \frac{\partial J}{\partial \boldsymbol{\Gamma}}
$$

Both terms are real. The explicit term is not a correction to be neglected: for
angle of attack it carries the freestream dependence of the Kutta–Joukowski
force **and** the body-to-stability-axis rotation, which together are a
first-order part of $dC_L/d\alpha$.

```text
dJ/dalpha = explicit fixed-circulation force and axis term
          + lambda^T db/dalpha
```

## `TotalDerivative`: the smallest entry point

```rust
use zephyr_vlm::gradients::{AeroOutput, TotalDerivative};

solver.solve(condition)?;
let lift = TotalDerivative::compute(solver, condition, AeroOutput::Lift)?;
println!("dCL/dalpha = {}", lift.by_alpha);
println!("dCL/dbeta  = {}", lift.by_beta);
```

`AeroOutput` is the scalar objective: `Lift`, `Drag`, `Side`, `Roll`, `Pitch`,
`Yaw`, or `InducedDrag` (the far-field Trefftz value). Objectives are this enum
rather than an integer selector, so a wrong objective is a compile error.

## `GradientCalculator`: the working interface

Install geometry once, choose what you want, then ask for a gradient set.

```rust
use zephyr_vlm::gradients::{AeroOutput, GradientCalculator, OutputSelection};

let mut gradients = GradientCalculator::new();
gradients.set_aircraft(aircraft);
gradients.set_output_selection(OutputSelection::only(AeroOutput::Lift));

let flight = gradients.compute_flight_condition_gradients(condition)?;
println!("dCL/dalpha = {}", flight.by_alpha.coefficients.lift);
println!("dCL/dM     = {}", flight.by_mach.coefficients.lift);
```

| Source | Call |
| --- | --- |
| `Aircraft` sections | `set_aircraft` |
| `MeshInput` direct mesh | `set_mesh` (fallible) |
| Moment reference alone | `set_reference_point` |

### Flight-condition gradients

`compute_flight_condition_gradients` returns a `FullGradients`, which carries an
`AeroGradients` per input:

| Field | Derivative with respect to |
| --- | --- |
| `by_alpha`, `by_beta` | the flight angles, per radian |
| `by_roll_rate`, `by_pitch_rate`, `by_yaw_rate` | the dimensionless body rates |
| `by_mach` | Mach number |
| `by_moment_reference[0..3]` | each moment-reference coordinate |

Each `AeroGradients` holds the six coefficient derivatives plus the far-field
induced-drag derivative, and `AeroGradients::get(output)` indexes it by
`AeroOutput`.

The moment-reference derivative is worth a note, because it is easy to get
wrong by hand. Moving the reference point does two independent things: it
changes every moment arm, and — because the point is also the centre the body
rates rotate about — it shifts the rotational onset velocity and therefore the
circulation. The second term vanishes when all three rates are zero, so at a
non-rotating condition only the moment outputs respond.

`compute_alpha_gradients` and `compute_beta_gradients` compute one input's set.
`flight_condition_gradients_from` reuses a solution a caller already has, which
is how you get several gradient sets from one operating point while paying for
the factorization once.

### Section geometry gradients

For a section model, the design variables are chord, twist, and the three
leading-edge coordinates of each section:

```rust
let geometry = gradients.compute_geometry_gradients(condition)?;
for section in &geometry.sections {
    println!(
        "{} section {}: dCL/dtwist = {:.6} per degree",
        section.surface_name,
        section.section_index,
        section.by_twist.coefficients.lift,
    );
}
let chord_sensitivities = geometry.chord_gradients(AeroOutput::Lift);
```

`GeometryVariable` is `Chord`, `Twist`, `LeadingEdgeX`, `LeadingEdgeY`, or
`LeadingEdgeZ`; `SectionGeometryGradients::get(variable)` indexes them.
Twist derivatives are **per degree**, matching the geometry convention.
`compute_section_gradient` does one section alone.

The whole set costs one primal solve. The adjoints are solved in a batch against
the cached transpose factorization, and each section direction is then
contracted against them — the aerodynamic system is never re-solved per
variable.

### Mesh-node gradients

For a direct mesh, the design variables are the node coordinates themselves:

```rust
let per_surface = gradients.compute_mesh_node_gradients(condition)?;
let wing = &per_surface[0];
let d_lift = wing.lift.node(0, 3).expect("node in range");
println!("dCL/dnode(0,3) = {d_lift:?}");
```

A `SurfaceMeshNodeGradients` carries one `MeshNodeGradients` per primary output
— `lift`, `drag`, `side`, `roll`, `pitch`, `yaw`, `induced_drag` — each a
row-major grid of `Vec3` cotangents indexed like `SurfaceMesh::nodes`.
`mesh_node_gradients_from` reuses an existing solution.

Mesh-node gradients hold the reference quantities fixed. If your reference area
or chord also depends on the design variables, that chain rule is yours.

### Control gradients

```rust
for control in gradients.compute_control_gradients(condition)? {
    println!(
        "dCL/d{} = {:.6} per radian",
        control.name, control.by_deflection.coefficients.lift,
    );
}
```

Every declared command is reported, including near-field and Trefftz induced
drag, honoring the `OutputSelection`. An omitted command is evaluated at zero
and appended. Rotations follow the flight condition's entry order, so keep that
order fixed.

Mesh-node and section gradients both work at nonzero deflection, and the mesh
gradients include the hinge direction changing with the mesh.

## Paying only for what you want

`OutputSelection` is a set of seven booleans with `all()`, `none()`, and
`only(output)` constructors and a `contains` test. An unselected output's
gradient array is left **exactly** zero, and its adjoint solve and the quadratic
velocity-kernel work behind it are skipped.

That is a real saving, not a cosmetic filter: selecting lift alone on a problem
with seven outputs removes six adjoint solves per gradient set.

For finer control still — individual `(output, input)` pairs, work counters, and
half-wing symmetry — use the request interface in
[Selective Analysis](./selective-analysis.md).

## Restrictions

Two restrictions apply to **geometry** adjoints, and both are errors rather than
silent approximations:

- **Explicit mirroring only.** A lattice using `SymmetryConfig` image planes is
  rejected. Duplicate a wing with `y_duplicate` instead; the image-plane
  geometry path is not represented in the reverse chain.
- **No section drag polars.** A `CDCL` polar makes a strip's viscous drag a
  function of its own lift coefficient, and hence of the strip area and the
  moment arm of its quarter chord. Those tangents do not exist yet, so a
  viscous mesh returns an error rather than an inviscid-only derivative.

Flight-condition and control derivatives carry the viscous chain and have
neither restriction.

Hinge placement and control span extent are configuration, not design
variables. There are no derivatives with respect to them.

## Lower-level kernels

Everything above is built from public pieces, so an objective the library does
not provide can be built from the same parts rather than by finite differencing
around it.

| Module item | Gives |
| --- | --- |
| `AdjointSolver` | transposed solves, single and batched, against the primal factorization |
| `horseshoe_velocity_jacobians` and the `horseshoe_velocity_vjp_*` functions | VJPs of the regularized Biot–Savart kernel with respect to the bound endpoints and the field point |
| `horseshoe_velocity_pg_derivative` | the Prandtl–Glauert derivative of the kernel |
| `aic_row_gradient`, `aic_column_gradient`, `aic_element_gradient` | AIC sensitivities |
| `velocity_row_gradient`, `velocity_column_gradient` | the same for the force-velocity operator |
| `force_circulation_gradient`, `force_alpha_gradient`, `force_beta_gradient`, `force_induced_velocity_gradient`, `force_geometry_jvp` | near-field force derivatives |
| `trefftz_circulation_gradients`, `trefftz_geometry_jvp` | far-field derivatives |
| `panel_normal_jvp`, `interpolated_twist_jvp` | section-to-lattice directions |

`AdjointSolver::solve_batched` is the one to reach for when a single output is
differentiated against many inputs: it amortizes the triangular solves over a
matrix of right-hand sides.

## Distributed loads

Coefficients are not the only differentiable output. For an arbitrary functional
of the **distributed** strip forces and moments — the aeroelastic case — see
[Coupled Interfaces](./coupled-interfaces.md) and `StripLoadAdjoint`.
