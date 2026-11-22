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

//! Trim derivatives checked against independent perturbed equilibrium solves.
#[path = "common/mesh_trim.rs"]
mod fixture;
use std::collections::BTreeMap;
use zephyr_vlm::core::VecX;
use zephyr_vlm::solver::VortexLatticeSolver;
use zephyr_vlm::trim::{
    AeroQuantity, ControlMixer, MeshTrimSolver, TrimVariable, linearize_aerodynamics,
};

#[test]
fn converged_states_and_all_output_partials_match_reconverged_differences() {
    let mesh = fixture::mesh();
    let condition = fixture::condition();
    let problem = fixture::problem();
    let external = [5.0, 0.0, -900.0, 0.0, 20.0, 0.0];
    let targets = fixture::targets(&mesh, &condition, &problem, &external);
    let solution = MeshTrimSolver
        .solve(&mesh, &condition, &problem, &targets, &external)
        .unwrap();
    fixture::close(solution.condition.alpha, condition.alpha);
    fixture::close(
        solution.condition.control_deflections[0].1,
        condition.control_deflections[0].1,
    );
    assert!(
        solution
            .diagnostics
            .residuals
            .iter()
            .all(|r| r.abs() < 1e-6)
    );
    let exact = solution
        .compute_partials(&mesh, &AeroQuantity::ALL, true)
        .unwrap();
    let adjoint = solution
        .compute_partials_adjoint(&mesh, &AeroQuantity::ALL, true)
        .unwrap();
    for (key, block) in &exact {
        for (a, b) in block.iter().zip(&adjoint[key]) {
            fixture::close(*a, *b);
        }
    }
    assert!(
        !exact
            .keys()
            .any(|(_, input)| input == "alpha" || input == "control_elevator")
    );
    let cases = [
        ("Wing", 17),
        ("Tail", 28),
        ("dynamic_pressure", 0),
        ("mach", 0),
        ("beta", 0),
        ("roll_rate", 0),
        ("pitch_rate", 0),
        ("yaw_rate", 0),
        ("control_flap", 0),
        ("control_tab", 0),
        ("moment_ref_x", 0),
        ("moment_ref_y", 0),
        ("moment_ref_z", 0),
        ("reference_area", 0),
        ("reference_chord", 0),
        ("reference_span", 0),
        ("trim_target_Fz", 0),
        ("trim_target_My", 0),
        ("external_Fz", 0),
        ("external_My", 0),
    ];
    for (name, index) in cases {
        let h = if name == "dynamic_pressure"
            || name.starts_with("trim_target")
            || name.starts_with("external_")
        {
            1e-3
        } else {
            1e-5
        };
        let mut results = Vec::new();
        for delta in [h, -h] {
            let mut m = mesh.clone();
            let mut c = condition.clone();
            let mut t = targets.clone();
            let mut e = external;
            match name {
                "Wing" | "Tail" => {
                    let s = usize::from(name == "Tail");
                    m.surfaces[s].nodes[index / 3][index % 3] += delta;
                }
                "dynamic_pressure" => c.dynamic_pressure += delta,
                "mach" => c.mach += delta,
                "beta" => c.beta += delta,
                "roll_rate" => c.roll_rate += delta,
                "pitch_rate" => c.pitch_rate += delta,
                "yaw_rate" => c.yaw_rate += delta,
                "moment_ref_x" => m.reference.moment_reference_point.x += delta,
                "moment_ref_y" => m.reference.moment_reference_point.y += delta,
                "moment_ref_z" => m.reference.moment_reference_point.z += delta,
                "reference_area" => m.reference.area += delta,
                "reference_chord" => m.reference.chord += delta,
                "reference_span" => m.reference.span += delta,
                "trim_target_Fz" => t[0] += delta,
                "trim_target_My" => t[1] += delta,
                "external_Fz" => e[2] += delta,
                "external_My" => e[4] += delta,
                _ => {
                    c.control_deflections
                        .iter_mut()
                        .find(|(n, _)| name == format!("control_{n}"))
                        .unwrap()
                        .1 += delta
                }
            }
            results.push(MeshTrimSolver.solve(&m, &c, &problem, &t, &e).unwrap());
        }
        for q in AeroQuantity::ALL {
            fixture::close(
                exact[&(q.name().into(), name.into())][index],
                (results[0].values[q as usize] - results[1].values[q as usize]) / (2.0 * h),
            );
        }
        for (i, variable) in problem.variables.iter().enumerate() {
            fixture::close(
                exact[&(variable.name.clone(), name.into())][index],
                (results[0].states[i] - results[1].states[i]) / (2.0 * h),
            );
        }
    }
    assert!(exact[&("CL".into(), "dynamic_pressure".into())][0].abs() > 1e-5);
    let mut stale = mesh.clone();
    stale.reference.chord += 0.01;
    assert!(solution.compute_partials(&stale, &[], false).is_err());
}

#[test]
fn unconverged_residual_state_jacobian_and_dimensional_loads_are_exact() {
    let mesh = fixture::mesh();
    let condition = fixture::condition();
    let mut problem = fixture::problem();
    problem
        .variables
        .push(TrimVariable::new("control_flap", 0.02));
    problem
        .balances
        .push(zephyr_vlm::trim::TrimBalance::new(AeroQuantity::Fx, 100.0));
    problem.validate(&mesh).unwrap();
    let z = VecX::from_vec(vec![0.04, -0.03, 0.01]);
    let c = problem.condition(&condition, &z).unwrap();
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&mesh).unwrap();
    solver.solve(&c).unwrap();
    let aero = linearize_aerodynamics(&mesh, &solver, &c, &AeroQuantity::ALL, false).unwrap();
    let state = problem.state_partials(&aero);
    for col in 0..z.len() {
        let mut values = Vec::new();
        for delta in [1e-6, -1e-6] {
            let mut zz = z.clone();
            zz[col] += delta;
            let cc = problem.condition(&condition, &zz).unwrap();
            solver.solve(&cc).unwrap();
            values.push(zephyr_vlm::trim::aerodynamic_values(&mesh, &solver, &cc));
        }
        for row in 0..12 {
            fixture::close(state[(row, col)], (values[0][row] - values[1][row]) / 2e-6);
        }
    }
    for name in ["reference_area", "reference_chord", "reference_span"] {
        for row in 6..12 {
            fixture::close(aero.partials[name][(row, 0)], 0.0);
        }
    }
}

#[test]
fn linear_mixer_reduces_redundant_controls_and_differentiates_realized_commands() {
    let mesh = fixture::mesh();
    let condition = fixture::condition();
    let mut problem = fixture::problem();
    problem.variables[1].name = "control_pitch".into();
    problem.mixer = Some(ControlMixer {
        commands: vec!["control_pitch".into()],
        surfaces: vec!["elevator".into(), "tab".into()],
        matrix: vec![vec![1.0], vec![-0.4]],
        offsets: vec![0.01, 0.02],
        bounds: BTreeMap::from([("elevator".into(), [-0.3, 0.3])]),
    });
    let true_condition = problem
        .condition(&condition, &VecX::from_vec(vec![0.075, 0.05]))
        .unwrap();
    let targets = fixture::targets(&mesh, &true_condition, &problem, &[0.0; 6]);
    let solution = MeshTrimSolver
        .solve(&mesh, &condition, &problem, &targets, &[0.0; 6])
        .unwrap();
    let exact = solution
        .compute_partials(&mesh, &[AeroQuantity::CD], true)
        .unwrap();
    fixture::close(solution.states[1], 0.05);
    for input in ["trim_target_Fz", "dynamic_pressure", "control_flap"] {
        fixture::close(
            exact[&("control_elevator".into(), input.into())][0],
            exact[&("control_pitch".into(), input.into())][0],
        );
        fixture::close(
            exact[&("control_tab".into(), input.into())][0],
            -0.4 * exact[&("control_pitch".into(), input.into())][0],
        );
    }
    let mut bad = problem.clone();
    bad.mixer
        .as_mut()
        .unwrap()
        .bounds
        .insert("elevator".into(), [-0.001, 0.001]);
    assert!(
        MeshTrimSolver
            .solve(&mesh, &condition, &bad, &targets, &[0.0; 6])
            .is_err()
    );
}

#[test]
fn singular_infeasible_and_invalid_trim_cases_fail_clearly() {
    let mesh = fixture::mesh();
    let condition = fixture::condition();
    let problem = fixture::problem();
    let targets = fixture::targets(&mesh, &condition, &problem, &[0.0; 6]);
    let mut singular = problem.clone();
    singular.variables[1].name = "beta".into();
    let mut symmetric = condition.clone();
    symmetric.beta = 0.0;
    symmetric.roll_rate = 0.0;
    symmetric.yaw_rate = 0.0;
    assert!(
        MeshTrimSolver
            .solve(&mesh, &symmetric, &singular, &targets, &[0.0; 6])
            .unwrap_err()
            .to_string()
            .contains("singular")
    );
    let mut bounds = problem.clone();
    bounds.variables[0].bounds = [-0.01, 0.01];
    bounds.variables[0].initial_guess = 0.0;
    bounds.variables[1].bounds = [-0.001, 0.001];
    assert!(
        MeshTrimSolver
            .solve(&mesh, &condition, &bounds, &targets, &[0.0; 6])
            .is_err()
    );
    let mut limit = problem.clone();
    limit.options.max_iterations = 1;
    assert!(
        MeshTrimSolver
            .solve(&mesh, &condition, &limit, &targets, &[0.0; 6])
            .unwrap_err()
            .to_string()
            .contains("iteration limit")
    );
    let mut invalid = problem.clone();
    invalid.variables[1].name = "alpha".into();
    assert!(invalid.validate(&mesh).is_err());
    invalid = problem.clone();
    invalid.balances[0].scale = 0.0;
    assert!(invalid.validate(&mesh).is_err());
    invalid = problem.clone();
    invalid.variables.pop();
    assert!(invalid.validate(&mesh).is_err());
    assert!(
        MeshTrimSolver
            .solve(&mesh, &condition, &problem, &[f64::NAN, 0.0], &[0.0; 6])
            .is_err()
    );
    assert!(
        MeshTrimSolver
            .solve(&mesh, &condition, &problem, &targets, &[f64::INFINITY; 6])
            .is_err()
    );
    invalid = problem.clone();
    invalid.variables[0].name = "unknown".into();
    assert!(invalid.validate(&mesh).is_err());
}

#[test]
fn external_nonlinear_structural_feedback_has_consistent_total_derivatives() {
    use zephyr_vlm::trim::MeshTrimSolution;
    let base_mesh = fixture::mesh();
    let base_condition = fixture::condition();
    let problem = fixture::problem();
    let targets = fixture::targets(&base_mesh, &base_condition, &problem, &[0.0; 6]);
    let compliance = 5e-4;
    let coupling = |q: f64| {
        let mut displacement = 0.0_f64;
        for _ in 0..12 {
            let mut mesh = base_mesh.clone();
            let mut condition = base_condition.clone();
            condition.dynamic_pressure = q;
            for surface in &mut mesh.surfaces {
                for node in &mut surface.nodes {
                    node.z += displacement * (node.y / 3.0).powi(2);
                }
            }
            let solution = MeshTrimSolver
                .solve(&mesh, &condition, &problem, &targets, &[0.0; 6])
                .unwrap();
            let partials = solution
                .compute_partials(&mesh, &[AeroQuantity::Fx, AeroQuantity::CD], true)
                .unwrap();
            let by_displacement = |output: &str| {
                mesh.surfaces
                    .iter()
                    .map(|s| {
                        s.nodes
                            .iter()
                            .enumerate()
                            .map(|(i, n)| {
                                partials[&(output.into(), s.name.clone())][3 * i + 2]
                                    * (n.y / 3.0).powi(2)
                            })
                            .sum::<f64>()
                    })
                    .sum::<f64>()
            };
            let residual = displacement + 0.2 * displacement.powi(3)
                - compliance * solution.values[AeroQuantity::Fx as usize];
            let jacobian = 1.0 + 0.6 * displacement.powi(2) - compliance * by_displacement("Fx");
            if residual.abs() < 1e-11 {
                let dz_dq =
                    compliance * partials[&("Fx".into(), "dynamic_pressure".into())][0] / jacobian;
                let total_cd = partials[&("CD".into(), "dynamic_pressure".into())][0]
                    + by_displacement("CD") * dz_dq;
                return (displacement, solution, dz_dq, total_cd);
            }
            displacement -= residual / jacobian;
        }
        panic!("external structural equilibrium did not converge");
    };
    let (w, solution, dw, dcd): (f64, MeshTrimSolution, f64, f64) =
        coupling(base_condition.dynamic_pressure);
    assert!(w.abs() > 1e-4);
    assert!(
        solution
            .diagnostics
            .residuals
            .iter()
            .all(|r| r.abs() < 1e-6)
    );
    let h = 0.1;
    let p = coupling(base_condition.dynamic_pressure + h);
    let m = coupling(base_condition.dynamic_pressure - h);
    fixture::close(dw, (p.0 - m.0) / (2.0 * h));
    fixture::close(
        dcd,
        (p.1.values[AeroQuantity::CD as usize] - m.1.values[AeroQuantity::CD as usize]) / (2.0 * h),
    );
}

#[test]
fn sideslip_trim_with_a_coefficient_target_has_analytical_state_sensitivities() {
    use zephyr_vlm::trim::{MeshTrimProblem, TrimBalance};
    let mesh = fixture::mesh();
    let c = fixture::condition();
    let problem = MeshTrimProblem {
        variables: vec![TrimVariable::new("beta", 0.0)],
        balances: vec![TrimBalance::new(AeroQuantity::CY, 1.0)],
        ..MeshTrimProblem::default()
    };
    let target = fixture::targets(&mesh, &c, &problem, &[0.0; 6]);
    let s = MeshTrimSolver
        .solve(&mesh, &c, &problem, &target, &[0.0; 6])
        .unwrap();
    fixture::close(s.states[0], c.beta);
    let exact = s
        .compute_partials(&mesh, &[AeroQuantity::CY, AeroQuantity::CD], false)
        .unwrap();
    assert!(
        !exact
            .keys()
            .any(|(_, input)| input.starts_with("external_") || input == "beta")
    );
    let h = 1e-6;
    let p = MeshTrimSolver
        .solve(&mesh, &c, &problem, &[target[0] + h], &[0.0; 6])
        .unwrap();
    let m = MeshTrimSolver
        .solve(&mesh, &c, &problem, &[target[0] - h], &[0.0; 6])
        .unwrap();
    fixture::close(
        exact[&("beta".into(), "trim_target_CY".into())][0],
        (p.states[0] - m.states[0]) / (2.0 * h),
    );
    fixture::close(
        exact[&("CD".into(), "trim_target_CY".into())][0],
        (p.values[1] - m.values[1]) / (2.0 * h),
    );
    fixture::close(exact[&("CY".into(), "trim_target_CY".into())][0], 1.0);
}

#[test]
fn a_supplied_solver_keeps_its_settings_through_the_trim() {
    let mesh = fixture::mesh();
    let condition = fixture::condition();
    let problem = fixture::problem();
    let external = [5.0, 0.0, -900.0, 0.0, 20.0, 0.0];
    let targets = fixture::targets(&mesh, &condition, &problem, &external);
    let plain = MeshTrimSolver
        .solve(&mesh, &condition, &problem, &targets, &external)
        .unwrap();
    // A default solver is exactly `solve`; mirror reuse changes no bit.
    for mirror_reuse in [false, true] {
        let mut solver = VortexLatticeSolver::new();
        solver.set_mirror_reuse(mirror_reuse);
        let supplied = MeshTrimSolver
            .solve_with_solver(solver, &mesh, &condition, &problem, &targets, &external)
            .unwrap();
        assert_eq!(supplied.solver.mirror_reuse(), mirror_reuse);
        assert_eq!(supplied.states, plain.states);
        assert_eq!(
            supplied.values.map(f64::to_bits),
            plain.values.map(f64::to_bits)
        );
        assert!(!supplied.solver.symmetric_factorization_applies());
    }
    // The fixture is mirrored with symmetric deflections, so the converged
    // solver, whose factorization the sensitivities reuse, uses the block form.
    let mut solver = VortexLatticeSolver::new();
    solver.set_symmetric_factorization(true);
    let symmetric = MeshTrimSolver
        .solve_with_solver(solver, &mesh, &condition, &problem, &targets, &external)
        .unwrap();
    assert!(symmetric.solver.symmetric_factorization_applies());
    for (a, b) in symmetric.values.iter().zip(&plain.values) {
        assert!((a - b).abs() <= 1e-12 * b.abs().max(1.0), "{a} vs {b}");
    }
    let exact = plain
        .compute_partials(&mesh, &AeroQuantity::ALL, true)
        .unwrap();
    let fast = symmetric
        .compute_partials(&mesh, &AeroQuantity::ALL, true)
        .unwrap();
    // Each partial rounds at the size of its Jacobian row; a balanced
    // quantity's derivatives vanish identically as cancelling untrimmed load
    // derivatives, the size of the largest dimensional entry in the column.
    let magnitude = |filter: &dyn Fn(&(String, String)) -> bool| {
        exact
            .iter()
            .filter(|(key, _)| filter(key))
            .flat_map(|(_, b)| b.iter())
            .fold(0.0_f64, |m, v| m.max(v.abs()))
    };
    for (key, block) in &exact {
        let mut scale = magnitude(&|k| k.0 == key.0);
        if problem.balances.iter().any(|b| b.quantity.name() == key.0) {
            scale = scale.max(magnitude(&|k| {
                k.1 == key.1 && AeroQuantity::ALL[6..].iter().any(|q| q.name() == k.0)
            }));
        }
        for (a, b) in fast[key].iter().zip(block.iter()) {
            assert!((a - b).abs() <= 1e-12 * scale, "{key:?}: {a} vs {b}");
        }
    }
}
