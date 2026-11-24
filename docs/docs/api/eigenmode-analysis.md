---
sidebar_position: 12
---

# EigenmodeAnalysis

`dynamics::EigenmodeAnalysis` assembles a linearized rigid-body system from
stability derivatives, mass properties, and a flight state, solves the
nonsymmetric complex eigenproblem, and classifies the modes.

The guide sheet is [Dynamic Modes](../guide/dynamic-modes.md), which carries the
rate-scaling caveat you should read before taking an absolute frequency.

```rust
use zephyr_vlm::dynamics::{EigenmodeAnalysis, FlightState, MassProperties};

let modes = EigenmodeAnalysis::new()
    .analyze_longitudinal(aircraft, mass, derivatives, flight)?;
println!("stable: {}", modes.is_stable());
```

## Systems

| Method | Size | States |
| --- | --- | --- |
| `analyze_longitudinal` | 4 | $u$, $\alpha$, $q$, $\theta$ |
| `analyze_lateral` | 4 | $\beta$, $p$, $r$, $\phi$ |
| `analyze_full` | 12 | $u$, $v$, $w$, $p$, $q$, $r$, $\phi$, $\theta$, $\psi$, $x$, $y$, $z$ |

`longitudinal_matrix`, `lateral_matrix`, and `full_matrix` return the matrix
without solving it. The 12-state system embeds the two 4-state blocks and adds
the heading and position kinematic rows; it is **block-decoupled by
construction**, with no aerodynamic or inertial cross-coupling between the
longitudinal and lateral blocks and no coupling through products of inertia. Add
that yourself if you need it.

An analysis is rejected with `ZephyrError::Eigensolver` when the mass or a
principal inertia is not positive, the derivatives are not marked valid, or the
velocity or density is not positive and finite.

## `MassProperties`

`mass`, `center_of_gravity`, `i_xx`, `i_yy`, `i_zz`, and the three products of
inertia — **stored positive**, with the tensor assembled as

$$
\mathbf{I} = \begin{pmatrix}
 I_{xx} & -I_{xy} & -I_{xz} \\
-I_{xy} &  I_{yy} & -I_{yz} \\
-I_{xz} & -I_{yz} &  I_{zz}
\end{pmatrix}
$$

`inertia_tensor()` and `set_inertia_tensor()` convert both ways, so a tensor
from elsewhere round-trips. `inverse_inertia_tensor()` is fallible.
`is_valid()` tests positivity. `MassFileData::total_properties()` produces one
from an AVL mass file with the parallel-axis transfer done.

## `FlightState`

Airspeed, `alpha`, `beta`, the normalized rates, the Euler attitude (`phi`,
`theta`, `psi`), `density`, and `gravity`. `dynamic_pressure()` is
$\tfrac12\rho V^2$ from its own fields and is what scales the aerodynamic terms;
`body_velocity()` is the body-axis velocity vector.

## `Eigenmode`

| Accessor | Value |
| --- | --- |
| `damping()` | $\mathrm{Re}(\lambda)$ |
| `frequency()` | $\lvert\mathrm{Im}(\lambda)\rvert$ |
| `period()` | $2\pi/\omega$, zero for a real eigenvalue |
| `time_constant()` | the **half-life**, $-\ln 2/\mathrm{Re}(\lambda)$, not the $1/e$ time |
| `damping_ratio()` | $-\mathrm{Re}(\lambda)/\lvert\lambda\rvert$ |
| `is_stable()`, `is_oscillatory()` | tests |

Each mode carries its `eigenvalue`, a normalized complex `eigenvector`, and a
`mode_type`.

## `ModeType`

`ShortPeriod`, `Phugoid`, `DutchRoll`, `RollSubsidence`, `Spiral`, plus
`Unknown` and `Other`; `name()` gives a display string.

Classification is heuristic with **absolute** thresholds in rad/s, chosen for
conventionally-sized aircraft. Longitudinally, an oscillatory mode above
1 rad/s with damping ratio above 0.1 is short period, below 1 rad/s is phugoid,
and anything else is `Other`. Laterally, any oscillatory mode is Dutch roll, a
real mode with damping below $-0.5$ is roll subsidence, and any other real mode
is spiral.

A small UAV or a very large aircraft will be misclassified. The eigenvalues are
always right; only the label is a guess.

## `EigenmodeResults`

`modes` (deterministically sorted), the exact `system_matrix` the eigenvalues
came from, `state_count`, and the reference velocity, span, and chord. Helpers:
`short_period()`, `phugoid()`, `dutch_roll()`, `modes_of_type(..)`, and
`is_stable()` — true only when every mode is stable.

[`EigenmodeAnalysis` rustdoc](https://docs.rs/zephyr-vlm/latest/zephyr_vlm/dynamics/struct.EigenmodeAnalysis.html)
