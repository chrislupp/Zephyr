---
sidebar_position: 8
---

# Dynamic Modes

`EigenmodeAnalysis` assembles a linearized rigid-body system from stability
derivatives, mass properties, and a flight state, then solves the nonsymmetric
complex eigenproblem and classifies the result into the conventional modes.

```rust
use zephyr_vlm::dynamics::{EigenmodeAnalysis, FlightState, MassProperties};

let mass = MassProperties {
    mass: 900.0,
    i_xx: 800.0,
    i_yy: 1_400.0,
    i_zz: 1_900.0,
    ..MassProperties::default()
};
let flight = FlightState { velocity: 50.0, ..FlightState::default() };

let modes = EigenmodeAnalysis::new()
    .analyze_longitudinal(aircraft, &mass, derivatives, &flight)?;

for mode in &modes.modes {
    println!(
        "{:<16} lambda = {:+.5}{:+.5}i   zeta = {:.3}   T = {:.3}",
        mode.mode_type.name(),
        mode.eigenvalue.re,
        mode.eigenvalue.im,
        mode.damping_ratio(),
        mode.period(),
    );
}
```

## The three systems

| Call | Size | States |
| --- | --- | --- |
| `analyze_longitudinal` | 4 | $u$, $\alpha$, $q$, $\theta$ |
| `analyze_lateral` | 4 | $\beta$, $p$, $r$, $\phi$ |
| `analyze_full` | 12 | $u$, $v$, $w$, $p$, $q$, $r$, $\phi$, $\theta$, $\psi$, $x$, $y$, $z$ |

The 12-state system embeds the two 4-state blocks at the corresponding indices,
adds the heading kinematic row, and adds the body-to-earth rotation rows for the
position states. It is block-decoupled by construction: there is no aerodynamic
or inertial cross-coupling between the longitudinal and lateral blocks, and
products of inertia do not couple them. If you need that coupling, the system
matrix is yours to modify.

`longitudinal_matrix`, `lateral_matrix`, and `full_matrix` return the matrix
without solving it. `EigenmodeResults::system_matrix` is the exact matrix the
reported eigenvalues came from, so a result can be checked independently.

## Rate derivative scaling

This is the one convention to get right before reading an absolute frequency.

The solver's body rates are nondimensionalized against a **unit freestream and
one geometry length unit**, and `StabilityDerivatives::by_roll_rate` and friends
are derivatives with respect to that quantity — see
[Conventions](../getting-started/conventions.md).

When forming the dimensional system, `EigenmodeAnalysis` multiplies each rate
derivative by the conventional half-span or half-chord factor: $b/2V$ for the
lateral rates, $c/2V$ for pitch rate. In other words, it reads the derivatives
as being with respect to $\hat p = pb/2V$, $\hat q = qc/2V$, $\hat r = rb/2V$.

Those two readings coincide only when the relevant reference length equals two
geometry length units. For any other model, scale the rate derivatives yourself
before handing over a `StabilityDerivatives`, or treat the resulting mode
frequencies as relative rather than absolute. Eigenvalue *signs* and therefore
stability conclusions are unaffected by a positive rescaling.

## Mass properties

```rust
let mass = MassProperties {
    mass: 900.0,
    i_xx: 800.0, i_yy: 1_400.0, i_zz: 1_900.0,
    i_xz: 40.0,
    ..MassProperties::default()
};
```

Inertias are about the centre of gravity. The products of inertia are stored
**positive**, so the tensor is assembled with the negative signs applied for
you:

$$
\mathbf{I} =
\begin{pmatrix}
 I_{xx} & -I_{xy} & -I_{xz} \\
-I_{xy} &  I_{yy} & -I_{yz} \\
-I_{xz} & -I_{yz} &  I_{zz}
\end{pmatrix}
$$

`inertia_tensor()` and `set_inertia_tensor()` convert both ways, so a tensor
from elsewhere round-trips without you having to remember the sign.
`inverse_inertia_tensor()` is fallible and reports a singular tensor rather than
producing garbage. The mass reader in
[AVL Input Files](./avl-files.md) produces these directly from a `.mass` file.

`FlightState` carries airspeed, the angles, the normalized rates, the Euler
attitude, density, and gravity. `dynamic_pressure()` is $\tfrac12\rho V^2$ from
its own fields, and is what scales the aerodynamic terms in the system matrix.

An analysis is rejected with `ZephyrError::Eigensolver` when the mass or a
principal inertia is not positive, the derivatives are not marked valid, or the
velocity or density is not positive and finite.

## Modes

Each `Eigenmode` carries an eigenvalue, a normalized complex eigenvector, and a
classification:

| Accessor | Value |
| --- | --- |
| `damping()` | $\mathrm{Re}(\lambda)$ |
| `frequency()` | $\lvert\mathrm{Im}(\lambda)\rvert$, rad/s |
| `period()` | $2\pi/\omega$, zero for a real eigenvalue |
| `time_constant()` | the **half-life**, $-\ln 2/\mathrm{Re}(\lambda)$ |
| `damping_ratio()` | $-\mathrm{Re}(\lambda)/\lvert\lambda\rvert$ |
| `is_stable()` | $\mathrm{Re}(\lambda) < 0$ |
| `is_oscillatory()` | nonzero imaginary part |

`time_constant` is the time to halve, not the $1/e$ time; it differs from the
usual definition by a factor of $\ln 2$.

Classification is heuristic and deliberately simple:

- **Longitudinal.** An oscillatory mode with frequency above 1 rad/s and damping
  ratio above 0.1 is `ShortPeriod`; an oscillatory mode below 1 rad/s is
  `Phugoid`; anything else is `Other`.
- **Lateral.** Any oscillatory mode is `DutchRoll`. A real mode with damping
  below $-0.5$ is `RollSubsidence`; any other real mode is `Spiral`.

Those thresholds are absolute, in rad/s, and they were chosen for
conventionally-sized aircraft. A small UAV or a very large aircraft will be
misclassified. The eigenvalues are always correct; only the label is a guess, so
read `eigenvalue` rather than relying on the name when the sizing is unusual.

Modes are sorted deterministically, so two runs of the same case list them in
the same order. `EigenmodeResults` offers `short_period()`, `phugoid()`,
`dutch_roll()`, `modes_of_type(..)`, and `is_stable()`, the last being true only
when every mode is stable.

Runnable: `cargo run --example stability_analysis`. See
[Stability Modes](../examples/stability-modes.md).
