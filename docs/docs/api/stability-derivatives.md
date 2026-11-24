---
sidebar_position: 9
---

# StabilityDerivativeCalculator

`aero::StabilityDerivativeCalculator` reports the derivatives of all six
coefficients with respect to angle of attack, sideslip, the three body rates,
and every named control command — plus the neutral point and static margin.

It costs no additional factorization: the angle and rate derivatives are linear
combinations of the six unit responses the solve already produced.

```rust
use zephyr_vlm::aero::StabilityDerivativeCalculator;

solver.solve(condition)?;
let d = StabilityDerivativeCalculator.compute(aircraft, solver, condition)?;
println!("CL_alpha = {:.5} /rad", d.by_alpha.lift);
println!("static margin = {:.4}", d.static_margin);
```

The solver must already hold a valid solution of exactly this condition;
otherwise the call returns `ZephyrError::InvalidSolution`.

## Methods

| Method | Returns |
| --- | --- |
| `compute` | `StabilityDerivatives`: stability-axis coefficient groups, controls, neutral point, static margin |
| `compute_body_axis` | `BodyAxisDerivatives`: body-axis force and moment **vectors** per input, without the stability-axis rotation and without the roll and yaw sign flips |

## `StabilityDerivatives`

Grouped by input, each group a full `Coefficients`:

| Field | With respect to |
| --- | --- |
| `by_alpha`, `by_beta` | the flight angles, per radian |
| `by_roll_rate`, `by_pitch_rate`, `by_yaw_rate` | the dimensionless body rates |
| `control_derivatives` | one `ControlDerivatives` per named command, per radian |

So `by_alpha.lift` is $C_{L\alpha}$ and `by_beta.roll` is $C_{l\beta}$.
Grouping by input rather than naming each derivative avoids the `CL`/`Cl`
collision, and makes a loop over inputs a loop over fields.

Also carried: `is_valid`, the reference quantities, `velocity`, the `alpha` and
`beta` the derivatives were taken at, `neutral_point_x`, and `static_margin`.

$$
x_{\text{np}} = x_{\text{ref}} - \frac{C_{m\alpha}}{C_{L\alpha}}c_{\text{ref}},
\qquad
\text{SM} = \frac{x_{\text{np}} - x_{\text{ref}}}{c_{\text{ref}}}
$$

Both are left at zero when $|C_{L\alpha}|$ is below the comparison epsilon.

## `ControlDerivatives`

`name`, the `deflection` it was evaluated at (summed over repeated entries),
`gain` (unity when unconfigured), `coefficients`, `body_force`, and
`hinge_moment` — which is **not computed** and is always zero.

Every command the lattice knows about is reported, whether or not it appears in
the flight condition's deflection list.

## Rate normalization

The rate derivatives are per unit of the solver's dimensionless rate, which is
an angular velocity against a unit freestream and one geometry length unit — not
the conventional $pb/2V$. `EigenmodeAnalysis` applies the conventional $b/2V$
and $c/2V$ factors when forming its dimensional system. Read the note in
[Dynamic Modes](../guide/dynamic-modes.md) before taking an absolute mode
frequency.

## When to use the gradient module instead

This is the conventional report. For an optimizer, `GradientCalculator` covers a
superset of the flight-condition inputs — Mach and the moment-reference
coordinates too — and respects an `OutputSelection` so you can pay for lift
alone. See [GradientCalculator](./gradient-calculator.md).

[`StabilityDerivativeCalculator` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/aero/struct.StabilityDerivativeCalculator.html)
