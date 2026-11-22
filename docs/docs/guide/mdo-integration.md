---
sidebar_position: 5
---

# MDO Integration

Zephyr provides a native Rust aerodynamic solver and analytical gradients.
The optional `philote` feature exposes two explicit gRPC disciplines:
`MeshDiscipline` maps section parameters to a structured mesh and reference
quantities, and `VlmDiscipline` maps a mesh and operating condition to aerodynamic
functions. The VLM discipline can solve trim internally while retaining its
explicit interface.

## Explicit servers

```bash
cargo run --features philote --example philote_mesh_server
cargo run --features philote --example philote_vlm_server
```

Configure surface names, dimensions, duplication, and controls before setup.
Per-surface lists accept JSON arrays or comma-separated strings. Controls use
material chordwise hinge coordinates, active span strips, and scalar radian
commands named `control_<name>`. The incoming mesh includes structural
deformation before the control normal rotations are applied.

Without trim, the VLM inputs are surface nodes, `reference_area`,
`reference_chord`, `reference_span`, `moment_ref_x/y/z`, `alpha`, `beta`, `mach`,
`dynamic_pressure`, normalized `roll_rate/pitch_rate/yaw_rate`, and each control
command. Its default outputs remain `CL,CD,CY,Cl,Cm,Cn`.

The selective `analysis_outputs`, `analysis_gradients` and `analysis_mode`
options described in the selective-analysis notes choose their own outputs and
partials, so they cannot be combined with `trim`, `outputs`, `fixed_inputs` or
`controls`.

The `outputs` option selects coefficients or dimensional resultants
`Fx,Fy,Fz,Mx,My,Mz`. Each is a scalar function with analytical gradients.
`fixed_inputs` suppresses selected input derivative columns while retaining
their physical effect. Geometry derivatives include moving control hinge axes.

## Solver accelerations

Two boolean options switch on the solver's opt-in accelerations for lattices
built from y-duplicated surfaces. Both default to `false`:

```json
{
  "surface_y_duplicated": [true, true],
  "mirror_reuse": true,
  "symmetric_factorization": true
}
```

- `mirror_reuse` derives the mirrored half of the influence operators by
  reflection instead of evaluating it. Results are bit-identical to having it
  off.
- `symmetric_factorization` factorizes two half-size systems instead of one
  full one when the system matrix is exactly mirror-symmetric. Outputs and
  partials then agree with the full factorization to rounding of the
  elimination, not bit for bit. Where the half-size form
  cannot apply, the discipline silently uses the full factorization and results
  are bit-identical to having it off: antisymmetric control deflections such as
  ailerons, surfaces that are not y-duplicated, and `analysis_mode` `half_wing`.

Each accepts a JSON boolean, `0`/`1`, or the strings `true`/`false` (also
`yes`/`no`); any other value is rejected. They apply to every solve the
discipline makes: ordinary function and partial evaluations, every iteration of
the internal trim and its sensitivities, and the selective `analysis_*` path.
They combine with every other option. Changing either clears the cached solved
state. Like `outputs`, an option omitted from a later request keeps its value.
Unlike the registry and trim options, they may also be changed after setup when
trim is configured; re-sending the current values there keeps the converged
trim. `VlmDiscipline::symmetric_factorization_applies` tells native callers
whether the cached solve used the half-size factorization.

## Internal trim configuration

The `trim` option accepts a JSON object or a JSON-encoded string. Set it in the
same request as the complete surface and control configuration, before setup:

```json
{
  "num_surfaces": 2,
  "surface_names": ["Wing", "Tail"],
  "surface_nx": [2, 2],
  "surface_ny": [3, 3],
  "surface_y_duplicated": [true, true],
  "controls": [
    {"surface": "Tail", "name": "elevator", "span_start": 0,
     "span_end": 3, "hinge": 1.2}
  ],
  "outputs": ["CL", "CD", "Fz", "My"],
  "trim": {
    "enabled": true,
    "variables": ["alpha", "control_elevator"],
    "balances": [
      {"quantity": "Fz", "scale": 1000.0},
      {"quantity": "My", "scale": 1000.0}
    ],
    "initial_guess": {"alpha": 0.03, "control_elevator": 0.0},
    "bounds": {"alpha": [-0.3, 0.3], "control_elevator": [-0.4, 0.4]}
  }
}
```

Selected `alpha`, `beta`, or control commands become outputs. Supply unselected
operating variables as inputs; supply no physical input for solved variables.
All solved variables are exposed even if `outputs` selects only one function.
Registry and trim settings are fixed after setup; create a new discipline to
change them.

Forces always act at bound-vortex centers, with every bound leg included in the
control-point normalwash. The `consistent_forces` option is still declared so
existing clients keep working: it accepts only `true`, which is the default and
changes nothing, and `false` is rejected because the legacy v0.2 force
convention was removed.

For the example, each evaluation includes scalar `trim_target_Fz` and
`trim_target_My`, plus all six `external_Fx/Fy/Fz/Mx/My/Mz` inputs. Each dimensional
balance is aerodynamic load plus prescribed external load minus target. For
example, external Fz of -1200 N and zero Fz target require an aerodynamic Fz of
1200 N. A custom `target_input` string can rename a balance's target variable.
Targets are differentiable physical inputs, so optimization can change them.

Coefficient balances such as `CL` and `Cm` use coefficient minus target and do
not add external loads. A problem using only coefficient balances declares no
external-load inputs. Mixed coefficient and dimensional rows are supported.

There must be one independent solved variable per unique balance quantity.
The example enforces Fz and My only; it does not establish all six aircraft
equilibrium equations. Speed, Mach, dynamic pressure, and rates are prescribed.
Atmosphere, propulsion, attitude transformations, and maneuver equations belong
to the surrounding model.

## Load conventions

Forces use mesh/body coordinate axes: x has positive drag and z has positive
lift at zero alpha. Moments are mathematical right-hand moments about
`moment_ref_x/y/z`. Dimensional Mx and Mz remove the conventional Cl and Cn sign
flips. All existing moment coefficients use reference chord for normalization.
Fx and Fz rotate near-field CL/CD into body axes and multiply by dynamic
pressure and reference area. No Trefftz drag substitution is made.

External forces and moments must use those same axes and moment reference.
The external model must differentiate any coordinate or moment-reference
transformation it applies. Normalized rate inputs also affect onset velocity
about the moment reference, and their analytical effects are included.

## Linked controls

A fixed linear mixer can reduce redundant physical controls to independent
trim commands. Add to the trim object, with physical controls already defined:

```json
{
  "mixer": {
    "commands": ["control_pitch"],
    "surfaces": ["elevator", "tab"],
    "matrix": [[1.0], [-0.4]],
    "offsets": [0.01, 0.02],
    "bounds": {"elevator": [-0.3, 0.3], "tab": [-0.2, 0.2]}
  }
}
```

Replace the solved `control_elevator` variable with `control_pitch` in this
configuration. Every mixer command must be a solved variable whose name differs
from physical control names. Rows are physical controls; columns are independent
commands. The realized commands are `offsets + matrix * commands`. Offsets
default to zero. Matrix weights and offsets are configuration rather than
design variables. Physical surface gains and duplication signs still apply.

Mixer-driven surface commands become outputs alongside the independent
commands, including their trim derivatives. They cease to be prescribed
inputs. Commands compose normal rotations in alphabetical physical-name order.

## Numerical solution and analytical gradients

`initial_guess` defaults to zero for each variable. `bounds` defaults to
unbounded; `variable_scales` defaults to one. Residual scales default to one in
each quantity's units. Optional `max_iterations` defaults to 30,
`max_backtracks` to 25, `absolute_tolerance` and `relative_tolerance` to 1e-10,
and `rank_tolerance` to 1e-12.

Scaled damped Newton uses analytical derivatives, backtracking, variable and
realized-control bounds, and a scaled Jacobian singular-value check. Each row
must satisfy absolute plus relative tolerance, where relative tolerance uses
that row's initial scaled residual. Small steps alone do not imply convergence.
Failures report operating states and unscaled physical balance residuals.
Configured initial guesses make repeated evaluations reproducible; continuation
warm starts are not currently provided.

For physical inputs p, solved states z, trim residual r, and an output f:

```text
r(p, z) = 0
dz/dp = -solve(dr/dz, dr/dp)
df/dp = partial(f)/partial(p) + partial(f)/partial(z) dz/dp
```

`compute_partials` returns the derivative of the complete trimmed explicit
mapping. These component partials include both circulation and trim response.
For example, CL can depend on dynamic pressure through dimensional trim even
though its direct dynamic-pressure derivative is zero. Direct and adjoint
methods differentiate the converged equations; Newton iterations and numerical
initial guesses are not differentiation inputs.

Function and derivative requests share an exact-input solved-state cache.
Derivative requests may precede function requests. Failed trim solves clear
the valid cache and return an error rather than approximate derivatives.
`VlmDiscipline::trim_diagnostics` gives native callers the cached convergence
residuals, iteration count, and scaled singular-value ratio.

## External coupling and examples

An external aeroelastic framework supplies deformed meshes and external loads,
receives a locally trimmed response, and converges structural and other feedback
connections using these explicit derivatives. External loads are fixed inside
each Zephyr solve. Loads depending on solved operating variables must be
connected through the external model and converged there. Angle of attack is
not an aircraft-attitude model for gravity.

The explicit mapping requires a nonsingular local trim Jacobian. Ineffective
controls, infeasible bounds, or branch changes can prevent a differentiable
local solution. Bound-active optimal control allocation is not provided.
Distributed loads are available through native strip-load APIs; transporting
them through this Philote discipline is a separate extension.

Run the native example without transport:

```bash
cargo run --example mesh_trim
```

For the explicit wire example, start a fresh VLM server and run this in another
terminal:

```bash
cargo run --features philote --example philote_trim_client
```

The client configures a deformed wing/tail trim, discovers solved variables,
requests functions and analytical gradients, and changes an external load.
The core `MeshTrimSolver` and `MeshTrimSolution::compute_partials` provide the
same mapping for native Rust optimizers without Philote dependencies.
