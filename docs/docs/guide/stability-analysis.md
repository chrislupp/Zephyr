---
sidebar_position: 6
---

# Stability and Control Derivatives

`StabilityDerivativeCalculator` reports the derivatives of all six coefficients
with respect to angle of attack, sideslip, the three body rates, and every named
control command. It needs a solved solver and costs no additional
factorization: the angle and rate derivatives are linear combinations of the six
unit responses the solve already produced.

```rust
use zephyr_vlm::aero::StabilityDerivativeCalculator;

solver.solve(condition)?;
let derivatives = StabilityDerivativeCalculator.compute(aircraft, solver, condition)?;

println!("CL_alpha = {:.5} /rad", derivatives.by_alpha.lift);
println!("Cm_alpha = {:.5} /rad", derivatives.by_alpha.pitch);
println!("Cn_beta  = {:.5} /rad", derivatives.by_beta.yaw);
println!("Cl_p     = {:.5}", derivatives.by_roll_rate.roll);
println!("neutral point x = {:.4}", derivatives.neutral_point_x);
println!("static margin   = {:.4}", derivatives.static_margin);
```

The solver must already hold a valid solution of exactly this condition;
`compute` checks and returns `ZephyrError::InvalidSolution` otherwise.

## Grouping

Derivatives are grouped by **input**, each group being a full `Coefficients`:

| Field | Derivative with respect to |
| --- | --- |
| `by_alpha` | angle of attack, per radian |
| `by_beta` | sideslip, per radian |
| `by_roll_rate`, `by_pitch_rate`, `by_yaw_rate` | the dimensionless body rates |
| `control_derivatives` | one `ControlDerivatives` per named command, per radian |

So `derivatives.by_alpha.lift` is $C_{L\alpha}$ and `derivatives.by_beta.roll`
is $C_{l\beta}$. Grouping by input rather than naming each derivative avoids the
`CL`/`Cl` collision that plagues flat naming schemes, and it means a loop over
inputs is a loop over struct fields rather than over a stringly-typed table.

`compute_body_axis` returns the same information as `BodyAxisDerivatives`:
force and moment **vectors** per input, in body axes, without the stability-axis
rotation and without the conventional roll and yaw sign flips.

## How the angle derivatives are formed

The unit responses are with respect to the freestream components, so the angle
derivatives come from the chain rule through

$$
\hat{\mathbf{V}}_\infty =
(\cos\alpha\cos\beta,\; \sin\beta,\; \sin\alpha\cos\beta)
$$

giving, for example, $\partial u/\partial\alpha = -\sin\alpha\cos\beta$ and
$\partial w/\partial\alpha = \cos\alpha\cos\beta$. The rate derivatives are the
rate unit responses directly.

This is exact for the circulation response and for the force response at the
reported operating point. It is a derivative of the linear potential-flow model,
not of any physical nonlinearity — see [Validity and Limits](./limitations.md).

## Neutral point and static margin

$$
x_{\text{np}} = x_{\text{ref}} - \frac{C_{m\alpha}}{C_{L\alpha}}\,c_{\text{ref}},
\qquad
\text{SM} = \frac{x_{\text{np}} - x_{\text{ref}}}{c_{\text{ref}}}
$$

Both are computed only when $|C_{L\alpha}|$ is above the comparison epsilon, and
are left at zero otherwise. With `x` aft, a neutral point behind the moment
reference point is a positive static margin.

## Control derivatives

One `ControlDerivatives` per command the lattice knows about, whether or not it
appears in the flight condition's deflection list:

| Field | Meaning |
| --- | --- |
| `name` | the command |
| `deflection` | the deflection it was evaluated at, summed over repeated entries |
| `gain` | linkage gain; unity when unconfigured |
| `coefficients` | coefficient derivatives per radian |
| `body_force` | body-axis force derivative |
| `hinge_moment` | **not computed**; always zero |

Each command's derivative comes from an analytical circulation derivative
(`VortexLatticeSolver::control_circulation_derivative`) through the same force
kernel, so it costs one triangular solve per command rather than a re-solve.

## The other route to the same numbers

`StabilityDerivativeCalculator` is the conventional report. For an optimizer you
usually want the gradient module instead:

| Want | Use |
| --- | --- |
| A conventional derivative table for a report or an eigenmode analysis | `StabilityDerivativeCalculator` |
| One output's derivative with respect to one input, cheaply | `analysis::Analysis` with a `GradientRequest` |
| Flight-condition derivatives including Mach and the moment reference point | `GradientCalculator::compute_flight_condition_gradients` |
| Derivatives with respect to geometry | `GradientCalculator`, section or mesh-node |

`FullGradients` from `GradientCalculator` covers alpha, beta, all three rates,
Mach, and the three moment-reference coordinates, which is a superset of the
flight-condition inputs `StabilityDerivatives` reports. It also respects an
`OutputSelection`, so you can pay for lift and pitching moment alone.

## Feeding the dynamics

`StabilityDerivatives` is the input `EigenmodeAnalysis` consumes, together with
`MassProperties` and a `FlightState`. Read [Dynamic Modes](./dynamic-modes.md)
before you do that — in particular the note on how the rate derivatives are
scaled into the dimensional system.

The runnable example is `cargo run --example stability_analysis`; see
[Stability Modes](../examples/stability-modes.md).
