---
sidebar_position: 8
---

# TrimSolver

`TrimConstraints` maps operating variables and ordered controls to fixed or
aerodynamic equations. Presets cover level flight, trimmed flight, steady
turn, and steady sideslip. `TrimSolver::solve` returns `TrimResult`;
nonconvergence is a status rather than an exception.
