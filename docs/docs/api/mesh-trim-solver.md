---
sidebar_position: 11
---

# MeshTrimSolver

`trim::MeshTrimSolver` converges a square trim system on a direct mesh and then
exposes the **derivatives of the converged mapping**. Trim becomes a
differentiable function of mesh nodes, targets, external loads, and prescribed
operating variables, rather than an iteration an optimizer has to see inside.

The guide sheet is [Trim](../guide/trim.md); the worked example is
[Mesh Trim](../examples/mesh-trim.md).

```rust
use zephyr_vlm::trim::{AeroQuantity, MeshTrimProblem, MeshTrimSolver, TrimBalance, TrimVariable};

let problem = MeshTrimProblem {
    variables: vec![
        TrimVariable::new("alpha", 0.03),
        TrimVariable::new("control_elevator", 0.0),
    ],
    balances: vec![
        TrimBalance::new(AeroQuantity::Fz, 1000.0),
        TrimBalance::new(AeroQuantity::My, 1000.0),
    ],
    ..MeshTrimProblem::default()
};
let external = [0.0, 0.0, -1200.0, 0.0, 0.0, 0.0];
let solution = MeshTrimSolver.solve(mesh, condition, &problem, &[0.0, 0.0], &external)?;
```

## `MeshTrimProblem`

| Field | Contents |
| --- | --- |
| `variables` | `TrimVariable`s in Jacobian-column order |
| `balances` | `TrimBalance`s in Jacobian-row order |
| `mixer` | an optional `ControlMixer` for redundant controls |
| `options` | `MeshTrimOptions` |

`validate(&mesh)` checks the configuration, physical control ownership, and that
the system is square, before anything is solved. `owns(name)` tests whether the
problem solves a given variable.

## `TrimVariable`

`name` is `alpha`, `beta`, or `control_<name>` — a physical control or a mixer
command. Plus a reproducible `initial_guess`, inclusive `bounds`, and a positive
`scale` used by the Newton solve and the rank check.

Solved variables must not also be supplied as prescribed operating values.

## `TrimBalance`

A `quantity`, a `target_input` naming the public input that carries this row's
target, and a fixed positive residual `scale` in that quantity's units.

`AeroQuantity` is `CL`, `CD`, `CY`, `Cl`, `Cm`, `Cn`, `Fx`, `Fy`, `Fz`, `Mx`,
`My`, `Mz`, with `ALL`, `from_name`, `name`, and `units`.

- **Coefficient** rows are *coefficient minus target* and need no external
  loads.
- **Dimensional** rows are *aerodynamic load plus prescribed external load minus
  target*. Body axes, mathematical right-hand moments about the configured
  reference point — so `Mx` and `Mz` carry **none** of the conventional `Cl`/`Cn`
  sign flips. External loads must match that convention and that reference
  point.

Mixed rows are supported. There must be exactly one independent variable per
unique balance quantity.

## `ControlMixer`

A fixed linear linkage, `delta = offsets + matrix * commands`. `commands` are
independent solved names in matrix-column order; `surfaces` are physical
control names in matrix-row order; `offsets` are per-surface radians;
`bounds` optionally limits the **realized** surface commands.

Every command must itself be a solved variable with a name distinct from any
physical control. The matrix and the offsets are configuration, not design
variables. Both the independent commands and the realized surface commands
become outputs with trim derivatives.

## `MeshTrimOptions`

Defaults: `max_iterations` 30, `max_backtracks` 25, `absolute_tolerance` and
`relative_tolerance` `1e-10`, `rank_tolerance` `1e-12`.

These are numerical controls, never physical design inputs. The relative
tolerance is referenced to each row's own initial scaled residual, and a small
step alone is never accepted as convergence.

## Solving

| Method | Effect |
| --- | --- |
| `solve(mesh, condition, problem, targets, external)` | converge the system |
| `solve_with_solver(..., solver)` | the same, starting from a preconfigured `VortexLatticeSolver` |

`solve_with_solver` is how `set_mirror_reuse` and
`set_symmetric_factorization` reach every iteration of the Newton solve and the
returned solution's cached factorization.

A failure is an **error**, reporting the operating states and the unscaled
physical residuals. It never returns an approximate derivative.

## `MeshTrimSolution`

| Field | Contents |
| --- | --- |
| `condition` | the solved `FlightCondition`, realized surface commands included |
| `solver` | the converged solver with its factorization |
| `states` | the solved independent variables |
| `values` | all twelve `AeroQuantity` values, coefficients then dimensional |
| `problem` | the configuration that was solved |
| `state_partials` | output derivatives with respect to the independent states |
| `jacobian` | the scaled residual Jacobian |
| `diagnostics` | `TrimDiagnostics`: `iterations`, **unscaled** physical `residuals`, `singular_value_ratio` |

A very small `singular_value_ratio` means the controls are nearly redundant for
these balances.

## Sensitivities

| Method | Direction |
| --- | --- |
| `compute_partials(mesh, quantities, nodes)` | direct; cheaper with few inputs |
| `compute_partials_adjoint(...)` | transpose solves; cheaper with many inputs and few functions |

Both return a `TrimPartials`, a `BTreeMap<(String, String), Vec<f64>>` from
`(output, input)` to a flattened block. Pass `nodes: false` to skip mesh adjoint
work when geometry is held fixed.

Every solved trim state is an output as well as the functions you request. Input
blocks include prescribed operating variables, balance targets, the six external
loads when a dimensional balance is present, dynamic pressure, reference
quantities, and the mesh nodes.

**Every block includes the trim response.** So $C_L$ can depend on dynamic
pressure through a dimensional trim even though its direct dynamic-pressure
derivative is zero. Newton iterations and initial guesses are not
differentiation inputs; only the converged equations are differentiated.

## Linearization

`trim::aerodynamic_values` and `trim::linearize_aerodynamics` are the layer
underneath: an `AeroLinearization` holding all twelve values and input Jacobian
blocks with circulation already eliminated, with operating variables held fixed.
Those are **not** trimmed totals. Set `nodes` false inside Newton iterations to
avoid mesh adjoint work.

[`MeshTrimSolver` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/trim/struct.MeshTrimSolver.html)
