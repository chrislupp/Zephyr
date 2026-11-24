---
sidebar_position: 13
---

# GradientCalculator

`gradients::GradientCalculator` is the working gradient interface: install
geometry once, choose outputs, then ask for flight-condition, section-geometry,
mesh-node, or control derivatives. Everything is analytical; there is no
finite-difference path.

The guide sheet is
[Analytical and Adjoint Gradients](../guide/automatic-differentiation.md).

```rust
use zephyr_vlm::gradients::{AeroOutput, GradientCalculator, OutputSelection};

let mut gradients = GradientCalculator::new();
gradients.set_aircraft(aircraft);
gradients.set_output_selection(OutputSelection::only(AeroOutput::Lift));
let flight = gradients.compute_flight_condition_gradients(condition)?;
println!("dCL/dalpha = {}", flight.by_alpha.coefficients.lift);
```

## Setup

| Method | Effect |
| --- | --- |
| `set_aircraft(&Aircraft)` | section model as the design source |
| `set_mesh(&MeshInput)` | direct mesh as the design source; fallible |
| `set_reference_point(Vec3)` | move the moment reference without other geometry changes |
| `set_output_selection(OutputSelection)` | which outputs to compute |

## Requests

| Method | Returns |
| --- | --- |
| `compute_flight_condition_gradients` | `FullGradients` |
| `flight_condition_gradients_from(solver, ...)` | the same, reusing a solution you already have |
| `compute_alpha_gradients`, `compute_beta_gradients` | one input's `AeroGradients` |
| `compute_geometry_gradients` | `GeometryGradients` for every section |
| `compute_section_gradient(condition, surface, section, variable)` | one `AeroGradients` |
| `compute_control_gradients` | `Vec<ControlGradient>`, per radian |
| `compute_mesh_node_gradients` | `Vec<SurfaceMeshNodeGradients>` |
| `mesh_node_gradients_from(solver, ...)` | the same, reusing a solution |

The `_from` variants require the solver to hold a valid solution of the given
condition on this calculator's geometry. Use them when several gradient sets
come from one operating point: the factorization is paid for once.

## `AeroOutput` and `OutputSelection`

`AeroOutput` is `Lift`, `Drag`, `Side`, `Roll`, `Pitch`, `Yaw`, or
`InducedDrag` — the last being the far-field Trefftz value.

`OutputSelection` is seven booleans with `all()`, `none()`, `only(output)`, and
`contains(output)`. An unselected output's array is left **exactly zero**, and
its adjoint solve and the quadratic velocity-kernel work behind it are skipped.
Selecting lift alone removes six adjoint solves per gradient set.

## `FullGradients`

| Field | With respect to |
| --- | --- |
| `by_alpha`, `by_beta` | the flight angles, per radian |
| `by_roll_rate`, `by_pitch_rate`, `by_yaw_rate` | the dimensionless body rates |
| `by_mach` | Mach number |
| `by_moment_reference[0..3]` | each moment-reference coordinate |

Each is an `AeroGradients`: the six coefficient derivatives in `coefficients`
plus `induced_drag`, indexable with `get(output)`.

Moving the reference point does two independent things — it changes every moment
arm, and, because the point is also the centre the body rates rotate about, it
shifts the rotational onset velocity and therefore the circulation. The second
term vanishes when all three rates are zero, so at a non-rotating condition only
the moment outputs respond.

## Section geometry

`GeometryVariable` is `Chord`, `Twist`, `LeadingEdgeX`, `LeadingEdgeY`, or
`LeadingEdgeZ`. Twist derivatives are **per degree**.

`GeometryGradients::sections` is surface-major. Each
`SectionGeometryGradients` carries `surface_name`, `section_index`, and
`by_chord`, `by_twist`, `by_leading_edge_x/y/z`, with `get(variable)` for
indexed access. `chord_gradients(output)` and `twist_gradients(output)` flatten
one output across all sections, and `design_variable_count()` reports the size.

The whole set costs one primal solve: the adjoints are batched against the
cached transpose factorization and each section direction is contracted against
them, so the aerodynamic system is never re-solved per variable.

## Mesh nodes

`SurfaceMeshNodeGradients` holds `surface_name` and one `MeshNodeGradients` per
primary output — `lift`, `drag`, `side`, `roll`, `pitch`, `yaw`, `induced_drag`.
Each is `nx`, `ny`, and a row-major `nodes: Vec<Vec3>` laid out like
`SurfaceMesh::nodes`, with a checked `node(i, j)`.

Reference quantities are held fixed. If your reference area or chord also
depends on the design variables, that chain rule is yours.

## Controls

`ControlGradient` is a `name` and a `by_deflection: AeroGradients`, per radian.
Every declared command is reported; an omitted one is evaluated at zero and
appended. Rotations follow the flight condition's entry order, so keep that
order fixed across an optimization.

## Restrictions

**Geometry** adjoints reject two configurations, as errors rather than silent
approximations:

- lattices using `SymmetryConfig` image planes — duplicate explicitly with
  `y_duplicate`;
- meshes carrying section drag polars, because the strip-area and quarter-chord
  tangents the viscous chain needs do not exist yet.

Flight-condition and control derivatives have neither restriction.

Hinge placement and control span extent are configuration; there are no
derivatives with respect to them.

## Lower-level pieces

| Item | Gives |
| --- | --- |
| `TotalDerivative::compute` | explicit plus implicit alpha and beta derivatives of one `AeroOutput` |
| `AdjointSolver` | transposed solves, single and `solve_batched`, against the primal factorization |
| `horseshoe_velocity_jacobians`, `horseshoe_velocity_vjp_*` | VJPs of the regularized Biot–Savart kernel |
| `horseshoe_velocity_pg_derivative` | the Prandtl–Glauert derivative |
| `aic_row_gradient`, `aic_column_gradient`, `aic_element_gradient` | AIC sensitivities |
| `velocity_row_gradient`, `velocity_column_gradient` | the same for the force-velocity operator |
| `force_circulation_gradient`, `force_alpha_gradient`, `force_beta_gradient`, `force_induced_velocity_gradient`, `force_geometry_jvp` | near-field force derivatives |
| `trefftz_circulation_gradients`, `trefftz_geometry_jvp` | far-field derivatives |
| `panel_normal_jvp`, `interpolated_twist_jvp` | section-to-lattice directions |

`solve_batched` is the one to reach for when a single output is differentiated
against many inputs.

[`GradientCalculator` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/gradients/struct.GradientCalculator.html)
