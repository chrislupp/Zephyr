---
sidebar_position: 4
---

# Automatic Differentiation

VLM provides analytical gradient computation for use in gradient-based optimization. Gradients are computed using the adjoint method, which efficiently computes derivatives of scalar outputs with respect to many inputs.

## Gradient Calculator

The `GradientCalculator` provides a high-level interface for computing derivatives of aerodynamic coefficients with respect to flight conditions and geometry parameters.

### Flight Condition Derivatives

```cpp
GradientCalculator grad_calc;
grad_calc.SetAircraft(aircraft);

FullGradients gradients = grad_calc.ComputeFlightConditionDerivatives();
```

The returned `FullGradients` struct contains `AeroGradients` for each perturbation variable ($\alpha$, $\beta$, $p$, $q$, $r$, $M$), where each `AeroGradients` holds:
- `dCL`, `dCD`, `dCY` — force coefficient derivatives
- `dCl`, `dCm`, `dCn` — moment coefficient derivatives
- `dCDi` — induced drag derivative

### Geometry Gradients

Derivatives with respect to wing section geometry (chord, twist, leading edge position):

```cpp
GeometryGradients geo_grads = grad_calc.ComputeGeometryGradients();

for (const auto& section_grad : geo_grads.section_gradients) {
    std::cout << "Surface: " << section_grad.surface_name
              << " Section: " << section_grad.section_index << std::endl;
    std::cout << "  dCL/dchord = " << section_grad.chord_derivs.dCL << std::endl;
    std::cout << "  dCL/dtwist = " << section_grad.twist_derivs.dCL << std::endl;
}
```

### Mesh Node Gradients

For direct mesh-level sensitivity, derivatives can be computed with respect to mesh node positions:

```cpp
grad_calc.SetMesh(solver.GetMesh());
auto mesh_grads = grad_calc.ComputeMeshNodeGradients();
```

## Adjoint Solver

The `AdjointSolver` provides the low-level adjoint solve for custom objective functions:

```cpp
AdjointSolver adjoint(solver);

// Get force gradients w.r.t. circulation
ForceGradients fg = ComputeForceGradients(
    solver.GetMesh(), solution, condition, aircraft);

// Solve adjoint system for CL objective
VecX adjoint_vector = adjoint.SolveAdjoint(fg.dCL_dGamma);
```

The adjoint vector can then be used to compute total derivatives with respect to any input parameter through the chain rule.
