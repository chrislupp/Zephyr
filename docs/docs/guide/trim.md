---
sidebar_position: 7
---

# Trim

There are two trim solvers, aimed at different problems.

| | `TrimSolver` | `MeshTrimSolver` |
| --- | --- | --- |
| Geometry | `Aircraft` sections | `MeshInput` direct mesh |
| Equations | coefficient targets on operating variables | coefficient targets **or** dimensional force and moment balances |
| External loads | none | prescribed `Fx`–`Mz` enter the residual |
| Result | `TrimResult`, nonconvergence as a status | `MeshTrimSolution`, nonconvergence as an error |
| Sensitivities | none | exact derivatives of the converged map |
| Use it for | a trimmed operating point in an analysis | trim inside an optimization or a coupled solve |

If you are trimming in order to differentiate, you want the second one.

## `TrimSolver`: constraints and presets

A trim problem is a `TrimConstraints`: for each of the five operating variables,
and for each named control, a `VariableConstraint` giving a kind and a target.

| `ConstraintType` | Equation |
| --- | --- |
| `Free` | omit this variable and its equation |
| `Fixed` | hold the variable at `target` |
| `Lift`, `Side` | match that force coefficient to `target` |
| `RollMoment`, `PitchMoment`, `YawMoment` | match that moment coefficient to `target` |

The system has to be square: every constraint that is not `Free` contributes one
equation, and every non-`Free` variable is one unknown.

```rust
use zephyr_vlm::trim::{TrimPresets, TrimSolver};

let result = TrimSolver::new().solve(
    aircraft,
    solver,
    &TrimPresets::level_flight(0.5),
    None,
);

println!("{:?} in {} iterations", result.status, result.iterations);
println!("alpha = {:.4} deg", result.state.alpha.to_degrees());
println!("residual = {:.3e}", result.residual_norm);
```

| Preset | Problem |
| --- | --- |
| `level_flight(CL)` | solve alpha for a lift target; sideslip and all rates fixed at zero |
| `level_trimmed_flight(CL, elevator)` | the same, plus the named control solved for zero pitching moment |
| `steady_sideslip(CL, beta)` | a lift target at a fixed, nonzero sideslip |
| `steady_turn(CL, bank)` | lift target with pitch and yaw rate left free. The bank argument is accepted for signature parity and is **not used**; the preset does not impose a turn kinematic relation. |

The Jacobian comes from the analytical stability derivatives, so a Newton step
is a step on exact derivatives of the reported forces, including the induced
velocity contribution — which is why it converges quadratically rather than
limping. Only the coefficients some equation drives and the variables some
constraint activates are computed.

Nonconvergence is a normal outcome, not an exception. `TrimResult::status` is a
`TrimStatus`:

| Status | Meaning |
| --- | --- |
| `Converged` | residual norm under `tolerance` |
| `NotStarted` | the solve never ran |
| `NoConstraints` | every constraint was `Free` |
| `AerodynamicFailure` | a solve inside the iteration failed |
| `SingularJacobian` | the trim Jacobian could not be factored |
| `BoundsExceeded` | a variable left its `TrimOptions` limit |
| `MaximumIterations` | the iteration limit was reached first |

`TrimResult` always carries the last state, aerodynamic results, derivatives,
and a human-readable `diagnostic`, converged or not, so a failure is
inspectable.

`TrimOptions` defaults: 20 iterations, a residual tolerance of `1e-6`, alpha and
beta limits of 1.5 rad, a rate limit of 5.0, and no relaxation. Pass an
`OperatingVariables` as the fourth argument to start from a specific guess
rather than zero.

## `MeshTrimSolver`: trim you can differentiate

This is the solver built for optimization and for coupled simulation. It
converges a square residual system on a prescribed mesh and then exposes the
**derivatives of the converged mapping** — so trim is a differentiable function
of mesh nodes, targets, external loads, and prescribed operating variables, not
an iteration an optimizer has to see inside.

A problem is a `MeshTrimProblem`: independent variables, balance equations, an
optional control mixer, and numerical options.

```rust
use zephyr_vlm::trim::{
    AeroQuantity, MeshTrimProblem, MeshTrimSolver, TrimBalance, TrimVariable,
};

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

let solution = MeshTrimSolver.solve(mesh, &condition, &problem, &[0.0, 0.0], &external)?;
```

### Variables

A `TrimVariable` is named `alpha`, `beta`, or `control_<name>`, carries a
reproducible initial guess, inclusive bounds, and a positive scale used by the
Newton solve and the rank check. Solved variables must not also be supplied as
prescribed operating values.

### Balances

A `TrimBalance` names an `AeroQuantity` and a positive residual scale in that
quantity's own units.

- **Coefficient balances** — `CL`, `CD`, `CY`, `Cl`, `Cm`, `Cn` — are
  *coefficient minus target*. They need no external loads.
- **Dimensional balances** — `Fx`, `Fy`, `Fz`, `Mx`, `My`, `Mz` — are
  *aerodynamic load plus prescribed external load minus target*. The axes are
  body axes and the moments are mathematical right-hand moments about the
  configured reference point, so `Mx` and `Mz` carry none of the conventional
  `Cl`/`Cn` sign flips. External loads must use the same axes and the same
  reference point.

Mixed coefficient and dimensional rows are fine. There must be exactly one
independent variable per unique balance quantity, and the system must be square.
`MeshTrimProblem::validate(&mesh)` checks that, plus physical control ownership,
before anything is solved.

### Control mixers

A `ControlMixer` is a fixed linear linkage, `delta = offsets + matrix * commands`,
reducing redundant physical surfaces to independent trim commands. Rows are
physical control names; columns are independent commands, each of which must
itself be a solved variable with a name distinct from any physical control.
Optional per-surface bounds limit the *realized* command.

The matrix and the offsets are configuration, not design variables. Both the
independent commands and the realized surface commands become outputs with trim
derivatives.

### Numerics

A scaled, damped Newton iteration on analytical derivatives, with backtracking,
variable and realized-control bounds, and a singular-value check on the scaled
Jacobian. Each row must satisfy an absolute plus a relative tolerance, where the
relative part is referenced to that row's own initial scaled residual — a small
step alone is never accepted as convergence.

Defaults: 30 iterations, 25 backtracks, `1e-10` absolute and relative
tolerances, `1e-12` rank tolerance. Initial guesses are configuration, which
makes repeated evaluations reproducible; there is no continuation warm start.

A failed solve returns an error reporting the operating states and the unscaled
physical residuals. It never returns an approximate derivative.

`solve_with_solver` takes a preconfigured `VortexLatticeSolver`, so
`set_mirror_reuse` and `set_symmetric_factorization` reach every iteration of
the Newton solve and the returned solution's cached factorization.

### Reading the solution

| Field | Contents |
| --- | --- |
| `condition` | the solved `FlightCondition`, including realized surface commands |
| `solver` | the converged solver, factorization included |
| `states` | the solved independent variables |
| `values` | all twelve `AeroQuantity` values, coefficients then dimensional |
| `state_partials` | output derivatives with respect to the independent states |
| `jacobian` | the scaled residual Jacobian |
| `diagnostics` | iterations, **unscaled** physical residuals, and the scaled singular-value ratio |

### Sensitivities

For physical inputs $p$, solved states $z$, residual $r$, and an output $f$:

$$
r(p, z) = 0,
\qquad
\frac{dz}{dp} = -\left(\frac{\partial r}{\partial z}\right)^{\!-1}\frac{\partial r}{\partial p},
\qquad
\frac{df}{dp} = \frac{\partial f}{\partial p} + \frac{\partial f}{\partial z}\frac{dz}{dp}
$$

```rust
let partials = solution.compute_partials(
    mesh,
    &[AeroQuantity::CL, AeroQuantity::CD],
    true, // include mesh-node blocks
)?;

let d_cd_d_q = partials[&("CD".into(), "dynamic_pressure".into())][0];
let d_alpha_d_fz = partials[&("alpha".into(), "external_Fz".into())][0];
```

`compute_partials` returns a `TrimPartials`, a map from `(output, input)` to a
flattened block. Every solved trim state is an output as well as the functions
you requested. Input blocks include the prescribed operating variables, balance
targets, the six external loads when a dimensional balance is present, dynamic
pressure, reference quantities, and — unless you pass `false` — the mesh nodes.

Every block includes the trim response. That has a consequence worth stating:
`CL` can depend on dynamic pressure through a dimensional trim, even though its
*direct* dynamic-pressure derivative is zero. Newton iterations and initial
guesses are not differentiation inputs; only the converged equations are
differentiated.

`compute_partials_adjoint` computes the same thing with transpose solves, which
is the cheaper direction when there are many inputs and few functions.

### What this does not do

The explicit mapping needs a nonsingular local trim Jacobian. Ineffective
controls, infeasible bounds, or a branch change can prevent a differentiable
local solution, and the solver reports that rather than guessing. Bound-active
optimal control allocation is not provided.

Prescribed quantities stay prescribed: speed, Mach, dynamic pressure, and rates.
Atmosphere, propulsion, attitude and gravity, and the full six-equation
equilibrium belong to the surrounding model — `alpha` here is an aerodynamic
angle, not an aircraft attitude.

Examples: `cargo run --example trim_flight` for the constraint solver, and
`cargo run --example mesh_trim` for the differentiable one. For trim inside a
gRPC discipline, see [MDO Integration](./mdo-integration.md).
