use pegasus_vlm::{
    aero::ForceCalculator,
    core::Vec3,
    geometry::{MeshInput, SurfaceMesh},
    gradients::{
        AeroGradients, COEFFICIENT_ORDER, MeshLinearization, MeshNodeDirection, OutputSelection,
    },
    solver::{FlightCondition, VortexLatticeSolver},
};

fn input(mirror: bool, multiple: bool) -> MeshInput {
    let mut input = MeshInput::default();
    for s in 0..if multiple { 2 } else { 1 } {
        let mut surface = SurfaceMesh::new(format!("surface{s}"), 2, 3);
        surface.y_duplicated = mirror;
        for i in 0..=2 {
            for j in 0..=3 {
                let y = j as f64 * 0.7 + 0.13;
                let chord = 1.3 - 0.13 * y;
                surface
                    .set_node(
                        i,
                        j,
                        Vec3::new(
                            0.2 * y + chord * i as f64 / 2.0 + s as f64 * 2.2,
                            y,
                            0.09 * y + 0.03 * i as f64 * y + s as f64 * 0.4,
                        ),
                    )
                    .unwrap();
            }
        }
        input.surfaces.push(surface);
    }
    input.reference.area = 4.0;
    input.reference.chord = 1.2;
    input.reference.span = 5.0;
    input.reference.moment_reference_point = Vec3::new(0.2, 0.1, -0.1);
    input
}
fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 0.08,
        beta: 0.02,
        mach: 0.2,
        roll_rate: 0.01,
        pitch_rate: -0.02,
        yaw_rate: 0.015,
        ..Default::default()
    }
}
fn solve(input: &MeshInput, c: &FlightCondition) -> VortexLatticeSolver {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    solver.solve(c).unwrap();
    solver
}
fn values(input: &MeshInput, c: &FlightCondition) -> AeroGradients {
    let s = solve(input, c);
    let f =
        ForceCalculator::new().compute_with_references(&input.reference, s.mesh(), s.solution(), c);
    AeroGradients {
        coefficients: f.coefficients,
        induced_drag: f.far_field_induced_drag,
    }
}
fn direction(input: &MeshInput, phase: f64) -> MeshNodeDirection {
    MeshNodeDirection {
        surfaces: input
            .surfaces
            .iter()
            .enumerate()
            .map(|(s, m)| {
                m.nodes
                    .iter()
                    .enumerate()
                    .map(|(i, _)| {
                        Vec3::from_fn(|k, _| ((i * 3 + k + s + 1) as f64 + phase).sin() * 0.1)
                    })
                    .collect()
            })
            .collect(),
    }
}
fn dot(a: &MeshNodeDirection, b: &MeshNodeDirection) -> f64 {
    a.surfaces
        .iter()
        .zip(&b.surfaces)
        .flat_map(|(a, b)| a.iter().zip(b))
        .map(|(a, b)| a.dot(b))
        .sum()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() < tol * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn forward_reverse_and_independent_resolves_agree_for_all_coefficients() {
    for (mirror, multiple, x_scale) in [
        (false, false, 1.0),
        (true, false, 1.0),
        (false, true, 1.0),
        (false, true, 12.0),
    ] {
        let mut input = input(mirror, multiple);
        for surface in &mut input.surfaces {
            for p in &mut surface.nodes {
                p.x *= x_scale;
            }
        }
        let c = condition();
        let solver = solve(&input, &c);
        let linear = MeshLinearization::new(&solver).unwrap();
        let directions = vec![direction(&input, 0.0), direction(&input, 0.7)];
        let forward = linear.jvp(&directions, OutputSelection::all()).unwrap();
        let mut seeds = Vec::new();
        for (i, o) in COEFFICIENT_ORDER.iter().enumerate() {
            let mut seed = AeroGradients::default();
            match o {
                pegasus_vlm::gradients::AeroOutput::Lift => seed.coefficients.lift = 1.0,
                pegasus_vlm::gradients::AeroOutput::Drag => seed.coefficients.drag = 1.0,
                pegasus_vlm::gradients::AeroOutput::Side => seed.coefficients.side = 1.0,
                pegasus_vlm::gradients::AeroOutput::Roll => seed.coefficients.roll = 1.0,
                pegasus_vlm::gradients::AeroOutput::Pitch => seed.coefficients.pitch = 1.0,
                pegasus_vlm::gradients::AeroOutput::Yaw => seed.coefficients.yaw = 1.0,
                pegasus_vlm::gradients::AeroOutput::InducedDrag => seed.induced_drag = 1.0,
            }
            seeds.push(seed);
            assert_eq!(seed.get(*o), 1.0, "{i}");
        }
        let reverse = linear.vjp(&seeds).unwrap();
        for (d, fwd) in directions.iter().zip(forward) {
            let mut plus = input.clone();
            let mut minus = input.clone();
            let h = 1e-5;
            for ((p, m), ds) in plus
                .surfaces
                .iter_mut()
                .zip(&mut minus.surfaces)
                .zip(&d.surfaces)
            {
                for ((p, m), d) in p.nodes.iter_mut().zip(&mut m.nodes).zip(ds) {
                    *p += h * d;
                    *m -= h * d;
                }
            }
            let vp = values(&plus, &c);
            let vm = values(&minus, &c);
            for (i, &o) in COEFFICIENT_ORDER.iter().enumerate() {
                close(fwd.get(o), dot(&reverse[i], d), 2e-12);
                close(fwd.get(o), (vp.get(o) - vm.get(o)) / (2.0 * h), 2e-7);
            }
        }
        // Arbitrary combined output seed equals the same linear combination.
        let mut mix = seeds[0];
        mix.coefficients.pitch = 0.3;
        mix.induced_drag = -0.8;
        let mixed = linear.vjp(&[mix]).unwrap();
        close(
            dot(&mixed[0], &directions[0]),
            dot(&reverse[0], &directions[0]) + 0.3 * dot(&reverse[4], &directions[0])
                - 0.8 * dot(&reverse[6], &directions[0]),
            1e-12,
        );
        assert!(linear.vjp(&[]).unwrap().is_empty());
        assert!(linear.jvp(&[], OutputSelection::all()).unwrap().is_empty());
        assert!(
            linear
                .jvp(&[MeshNodeDirection::default()], OutputSelection::all())
                .is_err()
        );
        let mut invalid = mix;
        invalid.induced_drag = f64::NAN;
        assert!(linear.vjp(&[invalid]).is_err());
        let zeros = linear.vjp(&[AeroGradients::default()]).unwrap();
        assert_eq!(dot(&zeros[0], &directions[0]), 0.0);
    }
}

#[test]
fn exact_scaling_translation_and_state_invalidation() {
    let mut input = input(false, false);
    input.reference.moment_reference_point = Vec3::zeros();
    let c = FlightCondition {
        alpha: 0.08,
        ..Default::default()
    };
    let mut solver = solve(&input, &c);
    let linear = MeshLinearization::new(&solver).unwrap();
    let scale = MeshNodeDirection {
        surfaces: input.surfaces.iter().map(|s| s.nodes.clone()).collect(),
    };
    let translation = MeshNodeDirection {
        surfaces: input
            .surfaces
            .iter()
            .map(|s| vec![Vec3::new(0.3, -0.2, 0.1); s.nodes.len()])
            .collect(),
    };
    let d = linear
        .jvp(&[scale, translation], OutputSelection::all())
        .unwrap();
    let v = values(&input, &c);
    for (i, o) in COEFFICIENT_ORDER.into_iter().enumerate() {
        close(
            d[0].get(o),
            v.get(o) * if (3..6).contains(&i) { 3.0 } else { 2.0 },
            1e-10,
        );
        if !(3..6).contains(&i) {
            close(d[1].get(o), 0.0, 1e-12);
        }
    }
    solver.set_mesh(&input).unwrap();
    assert!(MeshLinearization::new(&solver).is_err());
}

#[test]
fn regularized_kernel_includes_radius_chain_rule_and_exact_homogeneity() {
    use pegasus_vlm::gradients::horseshoe_velocity_jacobians;
    use pegasus_vlm::solver::horseshoe_vortex_velocity;
    let p = Vec3::new(0.7, 0.2, 0.03);
    let a = Vec3::new(0.0, -0.4, 0.1);
    let b = Vec3::new(0.1, 0.6, -0.03);
    let core = 0.08;
    let pg = 0.83;
    for bound in [false, true] {
        let jac = horseshoe_velocity_jacobians(&p, &a, &b, pg, core, bound);
        let h = 1e-6;
        let fd = (horseshoe_vortex_velocity(&p, &a, &b, pg, core + h, bound)
            - horseshoe_vortex_velocity(&p, &a, &b, pg, core - h, bound))
            / (2.0 * h);
        assert!((fd - jac.core_radius).norm() < 1e-8);
        let v = horseshoe_vortex_velocity(&p, &a, &b, pg, core, bound);
        let euler =
            jac.field_point * p + jac.bound_left * a + jac.bound_right * b + jac.core_radius * core;
        assert!((euler + v).norm() < 1e-12);
    }
}

#[test]
fn strip_load_tangents_include_moving_moment_points_and_match_adjoint() {
    use pegasus_vlm::gradients::{StripLoadAdjoint, StripLoadCotangent};
    let input = input(true, false);
    let c = condition();
    let solver = solve(&input, &c);
    let n = solver.mesh().strip_count();
    let points = vec![Vec3::new(0.3, 0.2, -0.1); n];
    let dp = vec![Vec3::new(0.07, -0.03, 0.01); n];
    let d = direction(&input, 0.3);
    let context = MeshLinearization::new(&solver).unwrap();
    let forward = context
        .strip_load_jvp(std::slice::from_ref(&d), &points, std::slice::from_ref(&dp))
        .unwrap();
    let seed = StripLoadCotangent {
        force: (0..n)
            .map(|i| Vec3::new(0.3, 0.1 * i as f64, -0.2))
            .collect(),
        moment: vec![Vec3::new(-0.1, 0.2, 0.3); n],
    };
    let adjoint = StripLoadAdjoint::new(&solver, &c, &points).unwrap();
    let nodes = adjoint
        .mesh_node_cotangents(&input, std::slice::from_ref(&seed))
        .unwrap();
    let point = adjoint
        .moment_point_cotangents(std::slice::from_ref(&seed))
        .unwrap();
    let product = |a: &StripLoadCotangent| {
        a.force
            .iter()
            .zip(&seed.force)
            .chain(a.moment.iter().zip(&seed.moment))
            .map(|(a, b)| a.dot(b))
            .sum::<f64>()
    };
    close(
        product(&forward[0]),
        dot(&nodes[0], &d) + point[0].iter().zip(&dp).map(|(a, b)| a.dot(b)).sum::<f64>(),
        2e-12,
    );
    let value = |step: f64| {
        let mut perturbed = input.clone();
        for (s, ds) in perturbed.surfaces.iter_mut().zip(&d.surfaces) {
            for (p, d) in s.nodes.iter_mut().zip(ds) {
                *p += step * d;
            }
        }
        let solver = solve(&perturbed, &c);
        let mut result = StripLoadCotangent::zeros(n);
        let fs = VortexLatticeSolver::freestream_direction(&c);
        for (i, vortex) in solver.mesh().vortices().iter().enumerate() {
            let s = vortex.strip_index;
            let state = solver.solution();
            let v = fs
                + Vec3::new(
                    state.induced_velocity[(i, 0)],
                    state.induced_velocity[(i, 1)],
                    state.induced_velocity[(i, 2)],
                );
            let f = 2.0
                * state.circulation[i]
                * v.cross(&(vortex.bound_vortex_right - vortex.bound_vortex_left));
            result.force[s] += f;
            result.moment[s] += (vortex.control_point - points[s] - step * dp[s]).cross(&f);
        }
        product(&result)
    };
    close(
        product(&forward[0]),
        (value(1e-5) - value(-1e-5)) / 2e-5,
        2e-7,
    );
    let zero = MeshNodeDirection {
        surfaces: input
            .surfaces
            .iter()
            .map(|s| vec![Vec3::zeros(); s.nodes.len()])
            .collect(),
    };
    let moving = context
        .strip_load_jvp(&[zero], &points, std::slice::from_ref(&dp))
        .unwrap();
    assert!(moving[0].force.iter().all(|v| v.norm() == 0.0));
    close(
        product(&moving[0]),
        point[0].iter().zip(&dp).map(|(a, b)| a.dot(b)).sum(),
        1e-13,
    );
    assert!(context.strip_load_jvp(&[d], &[], &[]).is_err());
}
