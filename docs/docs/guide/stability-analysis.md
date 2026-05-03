---
sidebar_position: 3
---

# Stability Analysis

VLM provides tools for computing stability derivatives and analyzing dynamic stability modes.

## Stability Derivatives

The `StabilityDerivativeCalculator` computes all aerodynamic derivatives by solving for unit perturbation responses:

```cpp
StabilityDerivativeCalculator deriv_calc;
StabilityDerivatives derivs = deriv_calc.ComputeDerivatives(
    aircraft, solver, condition);
```

### Available Derivatives

The `StabilityDerivatives` struct contains derivatives with respect to:

| Perturbation | Derivatives |
|-------------|-------------|
| Angle of attack ($\alpha$) | `CL_alpha`, `CD_alpha`, `Cm_alpha`, ... |
| Sideslip ($\beta$) | `CL_beta`, `CY_beta`, `Cl_beta`, `Cn_beta`, ... |
| Roll rate ($p$) | `CL_p`, `Cl_p`, `Cn_p`, ... |
| Pitch rate ($q$) | `CL_q`, `Cm_q`, ... |
| Yaw rate ($r$) | `CY_r`, `Cl_r`, `Cn_r`, ... |
| Control surfaces | `CL_d`, `Cm_d`, `Cl_d`, ... (per surface) |

### Static Stability

The neutral point and static margin are computed automatically:

```cpp
std::cout << "Neutral point: " << derivs.neutral_point_x << std::endl;
std::cout << "Static margin: " << derivs.static_margin << std::endl;
```

### Body Axis Derivatives

For body-axis derivatives (using $C_X$, $C_Y$, $C_Z$ instead of $C_L$, $C_D$):

```cpp
BodyAxisDerivatives body_derivs = deriv_calc.ComputeBodyAxisDerivatives(
    aircraft, solver, condition);
```

## Eigenmode Analysis

Dynamic stability is analyzed by constructing linearized equations of motion and computing eigenvalues.

### Mass Properties

```cpp
MassProperties mass;
mass.mass = 1000.0;    // kg
mass.Ixx = 1000.0;     // kg*m^2
mass.Iyy = 2000.0;
mass.Izz = 2500.0;
mass.Ixz = 100.0;
```

### Flight State

```cpp
FlightState flight;
flight.velocity = 50.0;   // m/s
flight.density = 1.225;    // kg/m^3
```

### Longitudinal Modes

```cpp
EigenmodeAnalysis analysis;
EigenmodeResults results = analysis.AnalyzeLongitudinal(
    aircraft, mass, derivs, flight);

for (const auto& mode : results.modes) {
    std::cout << GetModeTypeName(mode.type) << ": ";
    if (mode.IsOscillatory()) {
        std::cout << "freq=" << mode.GetFrequency() << " rad/s, "
                  << "damping=" << mode.GetDampingRatio();
    } else {
        std::cout << "time_const=" << mode.GetTimeConstant() << " s";
    }
    std::cout << (mode.IsStable() ? " [STABLE]" : " [UNSTABLE]")
              << std::endl;
}
```

The longitudinal modes include:
- **Short Period** — fast pitch oscillation
- **Phugoid** — slow speed/altitude exchange

### Lateral-Directional Modes

```cpp
EigenmodeResults lat_results = analysis.AnalyzeLateralDirectional(
    aircraft, mass, derivs, flight);
```

The lateral-directional modes include:
- **Dutch Roll** — coupled yaw-roll oscillation
- **Roll Subsidence** — pure roll damping
- **Spiral Mode** — slow divergence/convergence in heading
