---
sidebar_position: 10
---

# TrimSolver

`trim::TrimSolver` is the Newton driver for a section model: coefficient targets
on operating variables and controls. It reports a trimmed state and gives no
sensitivities — for trim inside an optimization or a coupled solve, see
[MeshTrimSolver](./mesh-trim-solver.md).

The guide sheet is [Trim](../guide/trim.md).

```rust
use zephyr_vlm::trim::{TrimPresets, TrimSolver, TrimStatus};

let result = TrimSolver::new().solve(
    aircraft, solver, &TrimPresets::level_flight(0.5), None,
);
if result.status != TrimStatus::Converged {
    eprintln!("{} ({:?})", result.diagnostic, result.status);
}
```

## `TrimConstraints`

A `VariableConstraint` — a `ConstraintType` and a target — for each of `alpha`,
`beta`, `roll_rate`, `pitch_rate`, `yaw_rate`, plus an ordered `controls` list of
`(name, constraint)` pairs.

| `ConstraintType` | Equation |
| --- | --- |
| `Free` | omit this variable and its equation |
| `Fixed` | hold the variable at `target` |
| `Lift`, `Side` | match that force coefficient |
| `RollMoment`, `PitchMoment`, `YawMoment` | match that moment coefficient |

The system must be square. `constraint_count()` and
`OperatingVariables::variable_count()` report both sides of that.

## `TrimPresets`

| Preset | Problem |
| --- | --- |
| `level_flight(CL)` | solve alpha for a lift target; sideslip and rates fixed at zero |
| `level_trimmed_flight(CL, elevator)` | the same, plus the named control solved for zero pitching moment |
| `steady_sideslip(CL, beta)` | a lift target at a fixed, nonzero sideslip |
| `steady_turn(CL, bank)` | lift target with pitch and yaw rate left free. The bank argument is accepted for signature parity and is **not used**; no turn kinematic relation is imposed. |

## `solve`

```text
solve(&Aircraft, &mut VortexLatticeSolver, &TrimConstraints, Option<&OperatingVariables>) -> TrimResult
```

The last argument is an initial guess; `None` starts from zero. The Jacobian
comes from the analytical stability derivatives — exact derivatives of the
reported forces, including the induced velocity — so the iteration converges
quadratically. Only the coefficients some equation drives and the variables some
constraint activates are computed.

## `TrimResult` and `TrimStatus`

Nonconvergence is a normal outcome, not an error. **Read `status` first**:
`TrimResult` carries its last state either way, so a failed trim read as a
converged one hands you a plausible alpha with a large residual behind it.

| `TrimStatus` | Meaning |
| --- | --- |
| `Converged` | residual norm under `tolerance` |
| `NotStarted` | the solve never ran |
| `NoConstraints` | every constraint was `Free` |
| `AerodynamicFailure` | a solve inside the iteration failed |
| `SingularJacobian` | the trim Jacobian could not be factored |
| `BoundsExceeded` | a variable left its `TrimOptions` limit |
| `MaximumIterations` | the iteration limit was reached first |

`TrimResult` also holds `converged`, `iterations`, `residual_norm`, `state`, the
full `aerodynamic_results`, the `derivatives` at that point, and a
human-readable `diagnostic`.

## `TrimOptions`

Defaults: 20 iterations, `tolerance` `1e-6`, `max_alpha` and `max_beta` 1.5 rad,
`max_rate` 5.0, `relaxation` 1.0, `verbose` false.

[`TrimSolver` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/trim/struct.TrimSolver.html)
