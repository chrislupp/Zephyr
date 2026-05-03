---
sidebar_position: 1
---

# VLM Theory

The Vortex Lattice Method is a numerical technique for computing aerodynamic loads on lifting surfaces in potential flow.

## Horseshoe Vortex Elements

VLM discretizes lifting surfaces into a grid of trapezoidal panels. Each panel carries a horseshoe vortex consisting of a bound vortex segment along the panel quarter-chord and two semi-infinite trailing vortices extending downstream.

The circulation strength $\Gamma_j$ of each horseshoe vortex is the unknown to be determined.

## Biot-Savart Law

The velocity induced by a vortex filament at a point $\mathbf{r}$ is computed using the Biot-Savart law:

$$
\mathbf{v} = \frac{\Gamma}{4\pi} \frac{(\mathbf{r}_1 \times \mathbf{r}_2)(r_1 + r_2)}{r_1 r_2 (r_1 r_2 + \mathbf{r}_1 \cdot \mathbf{r}_2)}
$$

where $\mathbf{r}_1$ and $\mathbf{r}_2$ are position vectors from the filament endpoints to the field point.

## Flow Tangency Condition

At each panel control point (located at the three-quarter chord), the normal component of the total velocity must vanish:

$$
\sum_{j=1}^{N} A_{ij} \Gamma_j = -\mathbf{V}_\infty \cdot \hat{\mathbf{n}}_i
$$

where $A_{ij}$ is the aerodynamic influence coefficient matrix and $\hat{\mathbf{n}}_i$ is the panel normal vector.

## Force Computation

### Near-Field: Kutta-Joukowski Theorem

Panel forces are computed using the Kutta-Joukowski theorem:

$$
\mathbf{F}_i = \rho \Gamma_i (\mathbf{V}_i \times \boldsymbol{\ell}_i)
$$

where $\mathbf{V}_i$ is the local velocity and $\boldsymbol{\ell}_i$ is the bound vortex vector.

### Far-Field: Trefftz Plane Analysis

The induced drag is also computed using Trefftz plane analysis, which evaluates the kinetic energy of the trailing vortex wake in a plane far downstream:

$$
C_{D_i} = \frac{1}{q_\infty S_{\text{ref}}} \sum_i \Gamma_i w_i \Delta y_i
$$

where $w_i$ is the downwash at the Trefftz plane.
