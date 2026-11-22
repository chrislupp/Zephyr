/*
    Project Zephyr

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-2198.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/

//! Opt-in block-symmetric factorization: it must engage exactly where the
//! system matrix is block-symmetric under the mirror pairing, fall back to the
//! full LU everywhere else, agree with the full LU to near machine precision
//! where it engages, and leave the default path untouched.

use zephyr_vlm::aero::ForceCalculator;
use zephyr_vlm::core::{MatX, Scalar, SymmetryType, Vec3, VecX};
use zephyr_vlm::geometry::{MeshControlSurface, MeshInput, ReferenceQuantities, SurfaceMesh};
use zephyr_vlm::gradients::{AdjointSolver, AeroGradients, GradientCalculator};
use zephyr_vlm::solver::{FlightCondition, SymmetryConfig, VortexLatticeSolver};

fn reference() -> ReferenceQuantities {
    ReferenceQuantities {
        area: 10.0,
        chord: 1.0,
        span: 10.0,
        moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
    }
}

/// Swept, dihedralled wing with a symmetric flap and an antisymmetric
/// aileron, so normals have three non-zero components and both deflection
/// symmetries are available.
fn wing(nx: usize, ny: usize, duplicated: bool) -> MeshInput {
    let mut surface = SurfaceMesh::new("wing", nx, ny);
    surface.y_duplicated = duplicated;
    for i in 0..=nx {
        let chordwise = 0.5 * (1.0 - (std::f64::consts::PI * i as Scalar / nx as Scalar).cos());
        for j in 0..=ny {
            let span = j as Scalar / ny as Scalar;
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        0.6 * span + chordwise * (1.0 - 0.3 * span),
                        0.1 + 4.9 * span,
                        0.35 * span + chordwise * 0.15 * span,
                    ),
                )
                .unwrap();
        }
    }
    let hinge = nx as Scalar * 0.7;
    surface.controls.push(MeshControlSurface::trailing_edge(
        "flap",
        0..ny / 2,
        vec![hinge; ny + 1],
    ));
    let mut aileron = MeshControlSurface::trailing_edge("aileron", ny / 2..ny, vec![hinge; ny + 1]);
    aileron.duplicate_sign = -1.0;
    surface.controls.push(aileron);
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

fn condition(
    alpha: Scalar,
    beta: Scalar,
    rates: [Scalar; 3],
    mach: Scalar,
    controls: &[(&str, Scalar)],
) -> FlightCondition {
    FlightCondition {
        alpha,
        beta,
        roll_rate: rates[0],
        pitch_rate: rates[1],
        yaw_rate: rates[2],
        mach,
        dynamic_pressure: 0.5 * 1.225 * 25.0 * 25.0,
        control_deflections: controls
            .iter()
            .map(|(name, value)| ((*name).to_string(), *value))
            .collect(),
    }
}

fn symmetric() -> FlightCondition {
    condition(0.08, 0.0, [0.0, 0.02, 0.0], 0.0, &[])
}

fn asymmetric() -> FlightCondition {
    condition(0.06, 0.04, [0.03, 0.02, -0.015], 0.2, &[])
}

fn solver(input: &MeshInput, symmetry: SymmetryConfig, block: bool) -> VortexLatticeSolver {
    let mut solver = VortexLatticeSolver::new();
    solver.set_symmetry(symmetry);
    solver.set_symmetric_factorization(block);
    solver.set_mesh(input).unwrap();
    solver
}

/// `max |a - b| / max |b|` over matching entries.
fn relative(a: &[Scalar], b: &[Scalar]) -> Scalar {
    assert_eq!(a.len(), b.len());
    let scale = b.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let error = a
        .iter()
        .zip(b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    if scale == 0.0 { error } else { error / scale }
}

/// Every solution vector the solver returns, flattened.
fn solution_values(solver: &VortexLatticeSolver) -> Vec<Scalar> {
    let solution = solver.solution();
    let mut values = Vec::new();
    for column in [
        &solution.circulation,
        &solution.circulation_per_u,
        &solution.circulation_per_v,
        &solution.circulation_per_w,
        &solution.circulation_per_p,
        &solution.circulation_per_q,
        &solution.circulation_per_r,
    ] {
        values.extend(column.iter().copied());
    }
    values.extend(solution.induced_velocity.iter().copied());
    values
}

fn bits(values: &[Scalar]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn gradient_values(gradients: &AeroGradients) -> [Scalar; 7] {
    let c = &gradients.coefficients;
    [
        c.lift,
        c.drag,
        c.side,
        c.roll,
        c.pitch,
        c.yaw,
        gradients.induced_drag,
    ]
}

fn force_values(solver: &VortexLatticeSolver, condition: &FlightCondition) -> Vec<Scalar> {
    let results = ForceCalculator::new().compute_with_references(
        solver.reference_quantities(),
        solver.mesh(),
        solver.solution(),
        condition,
    );
    let c = &results.coefficients;
    vec![
        c.lift,
        c.drag,
        c.side,
        c.roll,
        c.pitch,
        c.yaw,
        results.far_field_induced_drag,
        results.far_field_lift,
    ]
}

/// Rodrigues rotation, as the solver applies control deflections.
fn rotate(vector: Vec3, axis: &Vec3, angle: Scalar) -> Vec3 {
    let (sin, cos) = angle.sin_cos();
    vector * cos + axis.cross(&vector) * sin + axis * axis.dot(&vector) * (1.0 - cos)
}

/// Control-point normals after the condition's deflections.
fn deflected_normals(solver: &VortexLatticeSolver, condition: &FlightCondition) -> Vec<Vec3> {
    let mesh = solver.mesh();
    mesh.vortices()
        .iter()
        .zip(mesh.control_influences())
        .map(|(vortex, influences)| {
            let mut normal = vortex.normal_at_control_point;
            for (name, deflection) in &condition.control_deflections {
                for influence in influences.iter().filter(|item| item.name == *name) {
                    if influence.hinge_axis.norm_squared() > 1e-24 {
                        normal = rotate(
                            normal,
                            &influence.hinge_axis.normalize(),
                            deflection * influence.gain,
                        );
                    }
                }
            }
            normal
        })
        .collect()
}

/// `||A x - b||∞ / (||A||∞ ||x||∞)` for the system the solver actually
/// solved, including control-normal rotations.
fn relative_residual(solver: &VortexLatticeSolver, condition: &FlightCondition) -> Scalar {
    let matrix = solver.solution_aic();
    let x = &solver.solution().circulation;
    let freestream = VortexLatticeSolver::freestream_direction(condition);
    let omega = Vec3::new(
        condition.roll_rate,
        condition.pitch_rate,
        condition.yaw_rate,
    );
    let reference = *solver.reference_point();
    let normals = deflected_normals(solver, condition);
    let b = VecX::from_iterator(
        x.len(),
        solver
            .mesh()
            .vortices()
            .iter()
            .zip(&normals)
            .map(|(vortex, normal)| {
                if vortex.affected_by_freestream {
                    -(freestream + omega.cross(&(vortex.control_point - reference))).dot(normal)
                } else {
                    0.0
                }
            }),
    );
    relative_matrix_residual(matrix, x, &b)
}

fn relative_matrix_residual(matrix: &MatX, x: &VecX, b: &VecX) -> Scalar {
    let residual = matrix * x - b;
    let norm_a = (0..matrix.nrows())
        .map(|row| matrix.row(row).iter().map(|v| v.abs()).sum::<Scalar>())
        .fold(0.0, Scalar::max);
    residual.amax() / (norm_a * x.amax())
}

/// Accuracy bound on every solution value (circulations, unit responses,
/// induced velocities) and every force coefficient, relative to the largest
/// magnitude in the set. Measured differences are below 1e-15.
const SOLUTION_TOLERANCE: Scalar = 1e-13;
/// Accuracy bound on gradients, which pass through one more solve and a
/// longer chain of contractions. Measured differences are below 1e-14.
const GRADIENT_TOLERANCE: Scalar = 1e-12;
/// Bound on the normwise relative residual of either factorization. Both
/// measure around 2e-16 on these meshes.
const RESIDUAL_TOLERANCE: Scalar = 4e-15;

/// Solves `condition` with and without the block factorization and checks
/// the block form engaged and matches the full LU.
fn assert_block_matches_full(
    label: &str,
    input: &MeshInput,
    symmetry: SymmetryConfig,
    condition: &FlightCondition,
) {
    let mut full = solver(input, symmetry, false);
    let mut block = solver(input, symmetry, true);
    full.solve(condition).unwrap();
    block.solve(condition).unwrap();
    assert!(
        block.symmetric_factorization_applies(),
        "{label}: expected the block factorization to engage"
    );
    assert!(!full.symmetric_factorization_applies());

    let error = relative(&solution_values(&block), &solution_values(&full));
    assert!(
        error <= SOLUTION_TOLERANCE,
        "{label}: solution differs by {error:e}"
    );
    let error = relative(
        block.solution().circulation.as_slice(),
        full.solution().circulation.as_slice(),
    );
    assert!(
        error <= SOLUTION_TOLERANCE,
        "{label}: circulation differs by {error:e}"
    );
    let error = relative(
        &force_values(&block, condition),
        &force_values(&full, condition),
    );
    assert!(
        error <= SOLUTION_TOLERANCE,
        "{label}: forces differ by {error:e}"
    );

    let residual_full = relative_residual(&full, condition);
    let residual_block = relative_residual(&block, condition);
    assert!(
        residual_full <= RESIDUAL_TOLERANCE && residual_block <= RESIDUAL_TOLERANCE,
        "{label}: residuals full {residual_full:e}, block {residual_block:e}"
    );
}

/// Solves `condition` with and without the setting and checks the block form
/// declined, which must leave every output bit-identical.
fn assert_falls_back(
    label: &str,
    input: &MeshInput,
    symmetry: SymmetryConfig,
    condition: &FlightCondition,
) {
    let mut full = solver(input, symmetry, false);
    let mut block = solver(input, symmetry, true);
    full.solve(condition).unwrap();
    block.solve(condition).unwrap();
    assert!(block.symmetric_factorization());
    assert!(
        !block.symmetric_factorization_applies(),
        "{label}: expected a fallback to the full LU"
    );
    assert_eq!(
        bits(&solution_values(&block)),
        bits(&solution_values(&full)),
        "{label}: a fallback must reproduce the full LU exactly"
    );
}

#[test]
fn the_setting_is_off_by_default_and_toggling_it_restores_the_default_bits() {
    let input = wing(4, 8, true);
    let mut default_solver = VortexLatticeSolver::new();
    default_solver.set_mesh(&input).unwrap();
    assert!(!default_solver.symmetric_factorization());
    assert!(!default_solver.symmetric_factorization_applies());

    let mut toggled = VortexLatticeSolver::new();
    toggled.set_mesh(&input).unwrap();
    toggled.set_symmetric_factorization(true);
    toggled.solve(&asymmetric()).unwrap();
    assert!(toggled.symmetric_factorization_applies());
    toggled.set_symmetric_factorization(false);
    assert!(!toggled.symmetric_factorization_applies());

    for condition in [
        asymmetric(),
        condition(0.05, 0.02, [0.0; 3], 0.1, &[("flap", 0.15)]),
    ] {
        default_solver.solve(&condition).unwrap();
        toggled.solve(&condition).unwrap();
        assert!(!toggled.symmetric_factorization_applies());
        assert_eq!(
            bits(&solution_values(&toggled)),
            bits(&solution_values(&default_solver))
        );
        assert_eq!(
            bits(toggled.aic().as_slice()),
            bits(default_solver.aic().as_slice())
        );
        assert_eq!(
            bits(toggled.solution_aic().as_slice()),
            bits(default_solver.solution_aic().as_slice())
        );
    }
}

#[test]
fn undeflected_mirrored_lattices_use_the_block_form_at_any_flight_condition() {
    let none = SymmetryConfig::default();
    for (label, input) in [
        ("even span", wing(4, 8, true)),
        ("odd span", wing(5, 7, true)),
        ("single chordwise row", wing(1, 5, true)),
    ] {
        assert_block_matches_full(&format!("{label}, symmetric"), &input, none, &symmetric());
        assert_block_matches_full(&format!("{label}, asymmetric"), &input, none, &asymmetric());
    }
}

#[test]
fn symmetric_deflections_keep_the_block_form() {
    let none = SymmetryConfig::default();
    let input = wing(4, 8, true);
    assert_block_matches_full(
        "flap, symmetric condition",
        &input,
        none,
        &condition(0.05, 0.0, [0.0, 0.01, 0.0], 0.0, &[("flap", 0.2)]),
    );
    assert_block_matches_full(
        "flap, asymmetric condition",
        &input,
        none,
        &condition(0.06, 0.04, [0.03, 0.02, -0.015], 0.2, &[("flap", -0.1)]),
    );
    // Zero aileron command leaves the normals mirror-symmetric.
    assert_block_matches_full(
        "flap with zero aileron",
        &input,
        none,
        &condition(0.06, 0.0, [0.0; 3], 0.0, &[("flap", 0.1), ("aileron", 0.0)]),
    );
}

#[test]
fn ground_effect_images_keep_the_block_form_when_the_matrix_has_it() {
    // A z image is itself mirror-symmetric about y = 0. The decision rests on
    // the assembled matrix, which here passes the exact check.
    let ground = SymmetryConfig {
        z_symmetry: SymmetryType::Symmetric,
        z_symmetry_plane: -0.4,
        ..SymmetryConfig::default()
    };
    assert_block_matches_full("z image", &wing(4, 8, true), ground, &asymmetric());
}

#[test]
fn systems_without_the_block_structure_fall_back_to_the_full_lu() {
    let none = SymmetryConfig::default();
    let input = wing(4, 8, true);

    // Antisymmetric deflection: the mirrored normals rotate the other way.
    assert_falls_back(
        "aileron",
        &input,
        none,
        &condition(0.06, 0.0, [0.0; 3], 0.0, &[("aileron", 0.2)]),
    );
    assert_falls_back(
        "flap and aileron",
        &input,
        none,
        &condition(
            0.06,
            0.03,
            [0.01, 0.0, 0.0],
            0.1,
            &[("flap", 0.1), ("aileron", -0.1)],
        ),
    );

    // No duplicate, so no pairing at all.
    assert_falls_back("not duplicated", &wing(4, 8, false), none, &asymmetric());

    // Half-wing model: one side plus a y symmetry image, nothing to pair.
    let half = SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        ..SymmetryConfig::default()
    };
    assert_falls_back("half wing", &wing(4, 8, false), half, &symmetric());

    // Only one of two surfaces duplicated: some panels have no partner.
    let mut partial = wing(4, 8, true);
    let mut tail = SurfaceMesh::new("tail", 2, 3);
    for i in 0..=2 {
        for j in 0..=3 {
            tail.set_node(
                i,
                j,
                Vec3::new(4.0 + 0.3 * i as Scalar, -1.0 + 0.5 * j as Scalar, 0.4),
            )
            .unwrap();
        }
    }
    partial.surfaces.push(tail);
    assert_falls_back("partially duplicated", &partial, none, &asymmetric());

    // Mirrored about y = 0.05: `2*plane - y` rounds, so some mirror-image
    // entries differ in their last bits and the exact check rejects the
    // matrix. See `offset_mirror_planes_are_decided_by_the_matrix`.
    assert_falls_back(
        "offset mirror plane",
        &offset_wing(0.05, 0.37),
        none,
        &asymmetric(),
    );
}

/// Duplicated about `y = plane`, with linear node spacing so every value
/// involved comes from correctly rounded IEEE operations and whether the
/// exact check passes is portable.
fn offset_wing(plane: Scalar, step: Scalar) -> MeshInput {
    let mut surface = SurfaceMesh::new("wing", 4, 6);
    surface.y_duplicated = true;
    surface.y_duplicate_plane = plane;
    for i in 0..=4 {
        for j in 0..=6 {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        0.1 * j as Scalar + 0.25 * i as Scalar,
                        plane + step * j as Scalar,
                        0.05 * j as Scalar,
                    ),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: reference(),
    }
}

#[test]
fn offset_mirror_planes_are_decided_by_the_matrix() {
    // Mirroring about a plane other than y = 0 rounds `2*plane - y`. Whether
    // that rounding breaks the exact entry equality depends on the numbers,
    // so the decision is left to the matrix check rather than to the plane
    // metadata. Both outcomes occur, and both are correct.
    let none = SymmetryConfig::default();
    assert_block_matches_full(
        "offset plane, exactly symmetric",
        &offset_wing(0.35, 1.0 / 3.0),
        none,
        &asymmetric(),
    );
    assert_falls_back(
        "offset plane, rounded",
        &offset_wing(0.3, 0.7),
        none,
        &asymmetric(),
    );
}

#[test]
fn the_active_factorization_follows_the_deflection_state() {
    let input = wing(4, 8, true);
    let mut block = solver(&input, SymmetryConfig::default(), true);
    let mut full = solver(&input, SymmetryConfig::default(), false);
    for (condition, engaged) in [
        (asymmetric(), true),
        (
            condition(0.06, 0.0, [0.0; 3], 0.0, &[("aileron", 0.2)]),
            false,
        ),
        (condition(0.06, 0.0, [0.0; 3], 0.0, &[("flap", 0.2)]), true),
        (
            condition(0.06, 0.0, [0.0; 3], 0.0, &[("aileron", -0.1)]),
            false,
        ),
        (symmetric(), true),
    ] {
        block.solve(&condition).unwrap();
        full.solve(&condition).unwrap();
        assert_eq!(block.symmetric_factorization_applies(), engaged);
        let error = relative(&solution_values(&block), &solution_values(&full));
        assert!(error <= SOLUTION_TOLERANCE, "{error:e}");
    }
}

#[test]
fn gradients_through_the_block_form_match_the_full_lu() {
    let input = wing(4, 8, true);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();

    for condition in [
        asymmetric(),
        condition(0.06, 0.04, [0.03, 0.02, -0.015], 0.2, &[("flap", 0.15)]),
    ] {
        let mut full = solver(&input, SymmetryConfig::default(), false);
        let mut block = solver(&input, SymmetryConfig::default(), true);
        full.solve(&condition).unwrap();
        block.solve(&condition).unwrap();
        assert!(block.symmetric_factorization_applies());

        // Angle, rate and Mach derivatives.
        let flatten = |solver: &VortexLatticeSolver| -> Vec<Scalar> {
            let gradients = calculator
                .flight_condition_gradients_from(solver, &condition)
                .unwrap();
            [
                gradients.by_alpha,
                gradients.by_beta,
                gradients.by_roll_rate,
                gradients.by_pitch_rate,
                gradients.by_yaw_rate,
                gradients.by_mach,
            ]
            .iter()
            .flat_map(gradient_values)
            .collect()
        };
        let error = relative(&flatten(&block), &flatten(&full));
        assert!(
            error <= GRADIENT_TOLERANCE,
            "flight gradients differ by {error:e}"
        );

        // Mesh-node gradients run the transposed (adjoint) solves.
        let nodes = |solver: &VortexLatticeSolver| -> Vec<Scalar> {
            calculator
                .mesh_node_gradients_from(solver, &condition)
                .unwrap()
                .iter()
                .flat_map(|surface| {
                    [&surface.lift, &surface.drag, &surface.roll, &surface.pitch]
                        .into_iter()
                        .flat_map(|gradient| gradient.nodes.iter().flat_map(|v| [v.x, v.y, v.z]))
                        .collect::<Vec<_>>()
                })
                .collect()
        };
        let error = relative(&nodes(&block), &nodes(&full));
        assert!(
            error <= GRADIENT_TOLERANCE,
            "node gradients differ by {error:e}"
        );

        // Control derivatives, both the circulation solve and the residual.
        for name in ["flap", "aileron"] {
            let by_block = block
                .control_circulation_derivative(&condition, name)
                .unwrap();
            let by_full = full
                .control_circulation_derivative(&condition, name)
                .unwrap();
            let error = relative(by_block.as_slice(), by_full.as_slice());
            assert!(
                error <= GRADIENT_TOLERANCE,
                "{name} derivative differs by {error:e}"
            );
        }
    }
}

#[test]
fn transpose_solves_match_the_full_lu_and_satisfy_the_transposed_system() {
    let input = wing(4, 8, true);
    for condition in [
        asymmetric(),
        condition(0.06, 0.0, [0.0; 3], 0.0, &[("flap", 0.2)]),
    ] {
        let mut full = solver(&input, SymmetryConfig::default(), false);
        let mut block = solver(&input, SymmetryConfig::default(), true);
        full.solve(&condition).unwrap();
        block.solve(&condition).unwrap();
        assert!(block.symmetric_factorization_applies());

        let n = block.mesh().vortex_count();
        // Deliberately asymmetric right-hand sides.
        let rhs = MatX::from_fn(n, 3, |row, column| {
            ((row * 7 + column * 13) % 11) as Scalar - 5.0 + 0.1 * column as Scalar
        });
        let by_block = AdjointSolver::new(&block)
            .unwrap()
            .solve_batched(&rhs)
            .unwrap();
        let by_full = AdjointSolver::new(&full)
            .unwrap()
            .solve_batched(&rhs)
            .unwrap();
        let error = relative(by_block.as_slice(), by_full.as_slice());
        assert!(
            error <= SOLUTION_TOLERANCE,
            "batched transpose differs by {error:e}"
        );

        let column = rhs.column(1).into_owned();
        let single = AdjointSolver::new(&block).unwrap().solve(&column).unwrap();
        let error = relative(single.as_slice(), by_block.column(1).as_slice());
        assert!(
            error <= SOLUTION_TOLERANCE,
            "single transpose differs by {error:e}"
        );

        let transposed = block.solution_aic().transpose();
        let residual = relative_matrix_residual(&transposed, &single, &column);
        assert!(
            residual <= RESIDUAL_TOLERANCE,
            "transpose residual {residual:e}"
        );
    }
}

/// 256 panels: above the parallel assembly gates, and under `parallel` the
/// two half factorizations run concurrently.
#[test]
fn large_mirrored_lattices_match_the_full_lu() {
    let input = wing(8, 16, true);
    assert_block_matches_full(
        "256 panels",
        &input,
        SymmetryConfig::default(),
        &asymmetric(),
    );
    assert_block_matches_full(
        "256 panels, flap",
        &input,
        SymmetryConfig::default(),
        &condition(0.04, 0.02, [0.0, 0.01, 0.0], 0.3, &[("flap", 0.1)]),
    );
}

#[test]
fn mirror_reuse_and_the_block_form_compose() {
    // Mirror reuse is bit-identical to plain assembly, so the block form must
    // see the same matrix and produce the same bits with or without it.
    let input = wing(4, 8, true);
    for condition in [
        asymmetric(),
        condition(0.06, 0.0, [0.0; 3], 0.0, &[("flap", 0.2)]),
    ] {
        let mut plain = solver(&input, SymmetryConfig::default(), true);
        let mut mirrored = solver(&input, SymmetryConfig::default(), true);
        mirrored.set_mirror_reuse(true);
        plain.solve(&condition).unwrap();
        mirrored.solve(&condition).unwrap();
        assert!(mirrored.symmetric_factorization_applies());
        assert_eq!(
            bits(&solution_values(&mirrored)),
            bits(&solution_values(&plain))
        );
    }
}

#[test]
fn reused_deflected_solves_return_the_bits_of_a_fresh_solve() {
    // A repeat solve at unchanged deflections reuses the deflected
    // factorization and its unit responses. That reuse must be invisible.
    let input = wing(4, 8, true);
    for block in [false, true] {
        let first = condition(0.05, 0.01, [0.0, 0.01, 0.0], 0.1, &[("flap", 0.2)]);
        let second = condition(0.09, -0.02, [0.02, 0.0, 0.01], 0.1, &[("flap", 0.2)]);
        let mut reused = solver(&input, SymmetryConfig::default(), block);
        reused.solve(&first).unwrap();
        reused.solve(&first).unwrap();
        reused.solve(&second).unwrap();
        let mut fresh = solver(&input, SymmetryConfig::default(), block);
        fresh.solve(&second).unwrap();
        assert_eq!(
            bits(&solution_values(&reused)),
            bits(&solution_values(&fresh))
        );

        // Moving the reference point changes the rate unit responses, so it
        // must not reuse the cached ones.
        let point = Vec3::new(0.5, 0.1, -0.05);
        reused.set_reference_point(point);
        reused.solve(&second).unwrap();
        let mut moved = solver(&input, SymmetryConfig::default(), block);
        moved.set_reference_point(point);
        moved.solve(&second).unwrap();
        assert_eq!(
            bits(&solution_values(&reused)),
            bits(&solution_values(&moved))
        );
    }
}
