# Explicit Philote trim implementation plan

Extend the static aerodynamic `VlmDiscipline` with an optional internal trim
solve. The discipline remains an `ExplicitDiscipline`, served by the existing
`ExplicitServer`. Its declared functions of interest return analytical
derivatives through the converged trim solution. Trim is disabled by default.

This plan builds on the differentiable control surfaces in
`feature/control-surfaces`. The extension is implemented on
`feature/explicit-philote-trim`; the stages and validation below define its
acceptance criteria. Initial guesses are reproducible configuration; continuation
warm starts remain a follow-up capability.

The formulation follows the residual-based differentiation principle in
[A Novel Approach to Trimming Geometrically Nonlinearly Deformed Aircraft for
Gradient-Based Optimization Problems](https://doi.org/10.2514/6.2024-2592),
particularly Sections II.C and II.D. Here, trim residuals stay inside Pegasus,
and the external framework receives an explicit mapping and its derivatives.

## Scope and public behavior

The first implementation supports prescribed deformed meshes, trim variables
selected from `alpha`, `beta`, and named control commands, and either selected
dimensional force and moment balances or coefficient targets. Mach number,
dynamic pressure, and normalized body rates remain prescribed inputs. A
linear control mixer supports linked or redundant surfaces.

The incoming mesh represents structural deformation before control actuation.
The existing normal-rotation control model and its analytical hinge-direction
derivatives remain the aerodynamic model used during trim.

The initial functions of interest are the six existing coefficients, six
dimensional body force and moment components, and the solved trim variables.
Every exposed physical output receives complete derivatives through trim.
Other functions may be composed downstream if they depend only on exposed
outputs and inputs. A function that needs hidden aerodynamic state must be
implemented as an additional discipline output, including its derivative.

Structural equilibrium, load transfer, atmosphere, propulsion models, aircraft
attitude kinematics, and maneuver equations remain external. Optimal control
allocation, bound-active optimization, and solving for airspeed are follow-up
capabilities. Distributed loads over Philote are a separate integration step
described below.

Existing trim-disabled input names, output defaults, values, and derivatives
must remain compatible. The native legacy `TrimSolver` API remains available.

## Explicit input and output contract

Configuration is fixed before `setup` and variable discovery. Enabling trim
changes variable ownership; it cannot be toggled during evaluation without
rebuilding the registry and clearing the cache.

| Quantity | Contract with trim enabled |
| --- | --- |
| Surface node arrays and reference geometry | Physical inputs with analytical partials |
| Unselected flight variables and controls | Prescribed physical inputs |
| Selected trim variables | Scalar outputs using the corresponding operating-variable names |
| External force and moment | Physical inputs in the documented body frame and about the configured moment reference |
| Balance targets | Physical inputs, so their sensitivities are available to optimization |
| Coefficients and dimensional loads | Selected physical outputs evaluated at converged trim |
| Initial guesses, tolerances, scales, and limits | Numerical configuration options |
| Iteration counts and convergence diagnostics | Solver diagnostics, not differentiable physical outputs |

Solved states must always be discoverable outputs, even when the user requests
only a subset of aerodynamic functions. No name can simultaneously identify a
physical input and a solved output. Mixer-driven surface commands are derived
quantities rather than independently prescribed commands; requested realized
deflections can be exposed as outputs with their chain-rule derivatives.

Extend `fixed_inputs` validation to the new physical inputs. A selected trim
state cannot be a fixed input. Omitting external derivative columns must never
suppress the state derivatives needed to form the internal trim Jacobian.

Use a typed `TrimOptions` representation in the core and parse the Philote
option into it. The following is a proposed schema, not an existing API:

```json
{
  "trim": {
    "enabled": true,
    "variables": ["alpha", "control_elevator"],
    "balances": [
      {"quantity": "Fz", "target_input": "trim_target_Fz", "scale": 1000.0},
      {"quantity": "My", "target_input": "trim_target_My", "scale": 1000.0}
    ],
    "initial_guess": {"alpha": 0.05, "control_elevator": 0.0},
    "bounds": {
      "alpha": [-0.2, 0.3],
      "control_elevator": [-0.4, 0.4]
    },
    "max_iterations": 30,
    "absolute_tolerance": 1e-10,
    "relative_tolerance": 1e-10
  },
  "outputs": ["CL", "CD", "Cm", "Fz", "My"]
}
```

Angles and commands are radians. Balance scales have the units of their
associated quantities and must be positive. Tolerances apply to scaled
residuals. Target values and external loads are supplied as evaluation inputs,
not embedded in configuration. Resolve exact parser compatibility with the
current `outputs` option during implementation.

The example solves two selected balances. It does not enforce axial-force,
lateral-force, rolling-moment, or yawing-moment equilibrium. Documentation must
identify precisely which balances a configuration enforces.

## Internal residual and derivative formulation

Let `p` contain physical inputs, `z` the selected operating variables, and
`Gamma` the circulation. At each trial trim state, solve the VLM system:

```text
A(p, z) Gamma = b(p, z)
```

Eliminate circulation with the existing aerodynamic solve and analytical
linearization. Define trim residuals using the resulting loads:

```text
r(p, z) = D^-1 [S (L_aero(p, z) + L_external) - target]
```

Here `L` stacks dimensional body forces and moments, `S` selects the configured
components, and `D` contains fixed residual scales. Coefficient-target rows use
the selected aerodynamic coefficient minus its target and do not add
dimensional external loads. Reject incompatible quantity and unit combinations.

Dimensional aerodynamic loads must follow a single documented body-frame
convention. Audit the conversion from the existing force results: their body
load fields are normalized, moment signs need explicit treatment, and current
moment coefficients use the reference chord. Include dynamic-pressure,
reference-size, frame-rotation, and moment-reference derivatives in the
conversion. Do not substitute Trefftz drag into a near-field force balance.

All output and residual derivatives in the following equations already include
the circulation response, while holding trim states fixed:

```text
J = dr/dz
dz/dp = -solve(J, dr/dp)
df/dp = partial(f)/partial(p) + partial(f)/partial(z) dz/dp
```

For a scalar function of interest, the equivalent trim adjoint is:

```text
solve(J^T, lambda) = partial(f)/partial(z)^T
df/dp = partial(f)/partial(p) - lambda^T dr/dp
```

Use linear solves, never an explicit matrix inverse. These are derivatives of
the complete trimmed mapping and are returned as the explicit discipline's
partials. A coefficient can acquire a dynamic-pressure derivative through
dimensional trim even though its untrimmed dynamic-pressure partial is zero.

Start with reusable derivative blocks for requested functions and residual
rows, sharing aerodynamic adjoints and the small trim-Jacobian factorization.
Use the direct formulation for trim-state outputs and the adjoint formulation
where output count makes it advantageous. Do not compute mesh-coordinate
gradients at every Newton iteration: iteration requires state derivatives;
input gradients are needed at the converged solution. Never differentiate
through iteration counts, line-search decisions, or the initial guess.

## Control mixing and feasibility

Support an optional fixed linear map from independent commands `u` to named
surface commands:

```text
delta = delta_offset + B u
```

The initial mixer coefficients and offsets are configuration, not design
variables. Apply the mixer before the existing deterministic control-rotation
ordering and propagate its derivatives. Preserve surface gains and mirrored
control signs without applying them twice.

Require one independent trim variable per residual row. Validate unknown and
duplicate names, collisions with mesh variables, dimensions, units, bounds,
finite values, and contradictory command ownership during setup. Check
numerical Jacobian rank during solution; a square configuration is not a
guarantee of control authority.

Bounds limit the admissible trim solution and line-search steps. They do not
replace balance equations or introduce a control-allocation optimizer. A
nonzero residual at a bound is a failed trim solve. Report singularity,
stagnation, exceeded iteration limits, and inability to find an admissible
step with the current state and balance residuals. Detect actual surface-limit
violations after control mixing as well as independent-command limit violations.

## Solver and cache behavior

Implement damped Newton with state and residual scaling, a residual-based line
search, and explicit convergence checks on every selected balance. Define the
relative convergence reference from the initial scaled residual and retain an
absolute criterion for each row; report unscaled physical residuals as well.
Step stagnation alone is not convergence. Recompute the Jacobian at the
accepted converged state before differentiating.

The converged cache contains exact physical inputs, the configuration identity,
mesh, solved flight condition, aerodynamic solution, trim residuals, and the
trim Jacobian factorization. Cache derivative blocks lazily. Preserve existing
factorization reuse where geometry, Mach, and controls permit it.

`compute_partials` must work before `compute`, perform the same trim solve, and
return derivatives at the same accepted state. Failed solves must not populate
the valid-result cache. Option changes invalidate both solutions and gradients.
Serialize access consistently with the existing discipline cache.

Use configured initial guesses by default for reproducible evaluations. An
optional continuation warm start may reuse the last successful state, but it
is numerical state, not a physical input. Document that multiple trim roots
can make continuation branch-dependent. Gradient verification uses controlled
initialization and the same local solution branch.

## External aeroelastic coupling

For each input mesh and set of prescribed external loads, Pegasus returns a
locally trimmed aerodynamic response. The surrounding framework converges
structural deformation, load transfer, and other feedback connections using
these explicit derivatives.

External loads are held fixed during a local Pegasus trim solve. If weight
orientation, propulsion, or other loads depend on solved operating variables,
the external model must connect those dependencies and converge that feedback.
Angle of attack must not implicitly stand in for aircraft attitude in a
gravity model. All moments must refer to the same point, including derivatives
of any external moment transfer.

This architecture requires a nonsingular trim Jacobian at fixed discipline
inputs. A coupled aircraft equilibrium can exist even when that local
elimination is singular; such cases cannot be represented reliably by this
explicit mapping. Report the limitation through solver diagnostics rather
than returning regularized sensitivities as exact gradients.

The current Philote discipline exposes coefficients rather than distributed
structural loads. A production structural coupling will need an additional
output extension using `StripLoadAdjoint`, dimensional scaling, documented
force/moment locations, and their derivatives. Keep that extension separately
reviewable. A small coupling verification can use an exposed resultant load
and a scalar nonlinear structural response.

## Implementation stages

1. **Extract the shared aerodynamic evaluator.** Refactor mesh construction,
   flight-condition construction, force evaluation, and derivative access out
   of `src/philote/vlm.rs` into a small internal helper, provisionally
   `src/philote/aero_evaluation.rs`. Preserve the transport feature boundary.
   Add dimensional resultant conversion and selected output evaluation using
   core aerodynamic/gradient utilities. Acceptance: all existing Philote
   behavior passes unchanged and dimensional load derivatives agree with
   numerical checks at untrimmed conditions.

2. **Add a transport-independent trim model.** Add proposed modules
   `src/trim/problem.rs` and `src/trim/mesh_solver.rs`, exporting the required
   types through `src/trim/mod.rs`. Implement typed variables, balances,
   targets, mixer validation, residual evaluation, and the complete state
   Jacobian. Reuse the exact `GradientCalculator` control and flight-condition
   derivatives. Do not base the new implementation on the legacy trim
   Jacobian, which omits required cross-couplings. Acceptance: residual and
   Jacobian checks pass at arbitrary unconverged states with multiple controls.

3. **Implement nonlinear solution and trimmed sensitivities.** Add scaled
   damped Newton, admissible-step handling, convergence diagnostics, and
   converged-state direct/adjoint differentiation, provisionally in
   `src/trim/sensitivity.rs`. Reuse current mesh adjoints, including deformed
   hinge-axis derivatives. Acceptance: standalone mesh trim converges and
   re-solved finite differences verify solved-state and load sensitivities.

4. **Integrate with the explicit Philote discipline.** Extend option parsing,
   registry setup, selected output declarations, `fixed_inputs`, caching,
   `compute`, and `compute_partials` in `src/philote/vlm.rs`. Keep the current
   explicit service and default server behavior. Add trim configuration to
   `examples/philote_vlm_server.rs` documentation. Acceptance: clients discover
   the correct input/output ownership and receive complete trimmed partials
   through the explicit RPCs, including derivative-first requests.

5. **Add examples and complete validation.** Add a direct-mesh trim example
   and a Philote client workflow demonstrating target changes and total
   derivatives. Add a small external nonlinear coupling regression, with
   the structural model confined to test/example code. Update `README.md`,
   the MDO integration guide, and `CHANGELOG.md`. Acceptance: examples run,
   documented equations match the implemented contract, and repository checks
   pass.

Keep core trim and derivative code free of Philote, gRPC, and async-runtime
dependencies. The feature remains optional and disabled by default. No changes
to the legacy native trim algorithm are required for this extension.

## Validation matrix

Add focused core tests in `tests/mesh_trim.rs` and explicit-interface tests in
`tests/philote_trim.rs`, sharing fixtures where useful.

| Case | Required evidence |
| --- | --- |
| Trim disabled | Existing values, selected outputs, partials, and registry behavior remain compatible |
| Two-variable longitudinal trim | Every requested balance converges and solved states are exposed |
| Multiple interacting controls | Full state/control Jacobian, including cross-couplings, matches finite differences |
| Deformed mesh | Central differences with trim reconverged match analytical node gradients, including hinge motion |
| Targets and external loads | Target/external-load perturbations produce correct solved-state and output derivatives |
| Dimensional loads | Correct dynamic-pressure, reference-size, angle transformation, and moment-reference derivatives |
| Prescribed operating inputs | Mach, unselected controls, beta, and rates propagate correctly through trim where applicable |
| Control mixer | Linked commands, offsets, gains, and mirror signs give correct realized deflections and gradients |
| Direct and adjoint derivatives | Both formulations agree and match independently re-solved finite differences |
| Singular or infeasible case | Clear failure without cached success or fabricated derivatives |
| Cache and options | Derivative-first calls work; all relevant input/configuration changes invalidate results |
| External feedback | A nonlinear toy structural/load loop converges, and assembled totals match perturb-and-reconverge results |

Use finite differences for verification only. Check a range of perturbation
sizes, converge trim tighter than the derivative-check error budget, and
avoid branch changes or active limits in smooth-gradient tests. Check both
direct core calls and explicit Philote transport behavior.

Run the repository-required checks for the implementation:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo test --release --all-features
cargo test --no-default-features
RUSTDOCFLAGS="-D warnings" cargo doc --all-features --no-deps
cargo package
cargo llvm-cov --all-features --fail-under-lines 92
```

When the Docusaurus guide changes, also run `npm ci` and `npm run build` from
`docs/`. Measure primal and gradient work on a representative controlled mesh
to confirm that Newton iterations do not trigger unnecessary mesh-gradient
assembly and repeated derivative requests reuse the converged solution.

Completion requires a working explicit Philote trim option, analytical
partials for every declared function of interest, passing validation, and a
documented runnable example. Merely wrapping a trim iteration around existing
untrimmed partials does not meet the derivative contract.
