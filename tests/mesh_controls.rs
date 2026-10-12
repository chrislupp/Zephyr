//! Independent re-solve checks for externally deformed meshes with controls.
use pegasus_vlm::aero::{AerodynamicResults, ForceCalculator};
use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use pegasus_vlm::gradients::{
    AeroOutput, GradientCalculator, OutputSelection, StripLoadAdjoint, StripLoadCotangent,
};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

const OUTPUTS: [AeroOutput; 7] = [
    AeroOutput::Lift,
    AeroOutput::Drag,
    AeroOutput::Side,
    AeroOutput::Roll,
    AeroOutput::Pitch,
    AeroOutput::Yaw,
    AeroOutput::InducedDrag,
];

fn mesh(consistent: bool, mirrored: bool) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", 3, 4);
    surface.y_duplicated = mirrored;
    for i in 0..=3 {
        for j in 0..=4 {
            let eta = j as f64 / 4.0;
            let x = i as f64 / 3.0;
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(
                        0.2 * eta + (1.0 - 0.3 * eta) * x,
                        0.2 + 3.0 * eta,
                        0.15 * eta * eta + 0.12 * x * eta,
                    ),
                )
                .unwrap();
        }
    }
    surface.controls.push(MeshControlSurface::trailing_edge(
        "elevator",
        0..4,
        vec![1.4, 1.5, 1.6, 1.7, 1.8],
    ));
    let mut aileron = MeshControlSurface::trailing_edge("aileron", 1..4, vec![2.1; 5]);
    aileron.duplicate_sign = -1.0;
    aileron.gain = -0.7;
    surface.controls.push(aileron);
    let mut input = MeshInput {
        consistent_forces: consistent,
        ..MeshInput::default()
    };
    input.surfaces.push(surface);
    input.reference.area = 5.0;
    input.reference.chord = 0.8;
    input.reference.span = 6.0;
    input.reference.moment_reference_point = Vec3::new(0.3, 0.1, -0.05);
    input
}

fn flight() -> FlightCondition {
    FlightCondition {
        alpha: 0.08,
        beta: 0.02,
        mach: 0.3,
        roll_rate: 0.003,
        pitch_rate: -0.005,
        yaw_rate: 0.004,
        control_deflections: vec![("elevator".into(), 0.11), ("aileron".into(), -0.07)],
        ..FlightCondition::default()
    }
}

fn solve(
    input: &MeshInput,
    condition: &FlightCondition,
) -> (VortexLatticeSolver, AerodynamicResults) {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(input).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result = ForceCalculator::new().compute_with_references(
        &input.reference,
        solver.mesh(),
        &solution,
        condition,
    );
    (solver, result)
}

fn values(result: &AerodynamicResults) -> [f64; 7] {
    let c = result.coefficients;
    [
        c.lift,
        c.drag,
        c.side,
        c.roll,
        c.pitch,
        c.yaw,
        result.far_field_induced_drag,
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-7 + 2e-5 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn control_derivatives_include_induced_loads_and_match_resolves() {
    for consistent in [false, true] {
        let input = mesh(consistent, true);
        let mut calculator = GradientCalculator::new();
        calculator.set_mesh(&input).unwrap();
        for condition in [
            flight(),
            FlightCondition {
                control_deflections: vec![],
                ..flight()
            },
            FlightCondition {
                control_deflections: vec![
                    ("elevator".into(), 0.05),
                    ("aileron".into(), 0.03),
                    ("elevator".into(), 0.04),
                ],
                ..flight()
            },
        ] {
            let gradients = calculator.compute_control_gradients(&condition).unwrap();
            assert_eq!(gradients.len(), 2);
            for gradient in gradients {
                let mut plus = condition.clone();
                let mut minus = condition.clone();
                if !plus
                    .control_deflections
                    .iter()
                    .any(|(n, _)| n == &gradient.name)
                {
                    plus.control_deflections.push((gradient.name.clone(), 0.0));
                    minus.control_deflections.push((gradient.name.clone(), 0.0));
                }
                for (n, v) in &mut plus.control_deflections {
                    if n == &gradient.name {
                        *v += 1e-6;
                    }
                }
                for (n, v) in &mut minus.control_deflections {
                    if n == &gradient.name {
                        *v -= 1e-6;
                    }
                }
                let p = values(&solve(&input, &plus).1);
                let m = values(&solve(&input, &minus).1);
                for (k, output) in OUTPUTS.iter().enumerate() {
                    close(gradient.by_deflection.get(*output), (p[k] - m[k]) / 2e-6);
                }
            }
        }
    }
}

#[test]
fn deformed_hinge_mesh_gradients_match_every_node_resolve() {
    for consistent in [false, true] {
        let input = mesh(consistent, true);
        let condition = flight();
        let mut calculator = GradientCalculator::new();
        calculator.set_mesh(&input).unwrap();
        let generic = calculator.compute_mesh_node_gradients(&condition).unwrap();
        calculator.set_output_selection(OutputSelection {
            induced_drag: false,
            ..OutputSelection::all()
        });
        let batched = calculator.compute_mesh_node_gradients(&condition).unwrap();
        for node in 0..input.surfaces[0].nodes.len() {
            for d in 0..3 {
                let mut plus = input.clone();
                let mut minus = input.clone();
                plus.surfaces[0].nodes[node][d] += 1e-6;
                minus.surfaces[0].nodes[node][d] -= 1e-6;
                let p = values(&solve(&plus, &condition).1);
                let m = values(&solve(&minus, &condition).1);
                for g in [&generic[0], &batched[0]] {
                    let blocks = [&g.lift, &g.drag, &g.side, &g.roll, &g.pitch, &g.yaw];
                    for (k, block) in blocks.iter().enumerate() {
                        close(block.nodes[node][d], (p[k] - m[k]) / 2e-6);
                    }
                }
                close(generic[0].induced_drag.nodes[node][d], (p[6] - m[6]) / 2e-6);
            }
        }
    }
}

fn functional(result: &AerodynamicResults, points: &[Vec3], column: &StripLoadCotangent) -> f64 {
    result
        .surface_loads
        .iter()
        .flat_map(|s| &s.strip_loads)
        .map(|s| {
            let force = s.area * s.body_force;
            let moment =
                s.area * s.chord * Vec3::new(-s.body_moment.x, s.body_moment.y, -s.body_moment.z)
                    + (s.quarter_chord_position - points[s.strip_index]).cross(&force);
            column.force[s.strip_index].dot(&force) + column.moment[s.strip_index].dot(&moment)
        })
        .sum()
}

#[test]
fn strip_load_pullbacks_include_controls_deformation_and_flight_angles() {
    let input = mesh(true, true);
    let condition = flight();
    let (solver, _) = solve(&input, &condition);
    let n = solver.mesh().strip_count();
    let points: Vec<_> = (0..n)
        .map(|j| Vec3::new(0.2, j as f64 * 0.4, 0.1))
        .collect();
    let columns: Vec<_> = (0..2)
        .map(|k| StripLoadCotangent {
            force: (0..n)
                .map(|j| Vec3::new(0.1 * j as f64, 0.3, 0.8 + k as f64))
                .collect(),
            moment: (0..n)
                .map(|j| Vec3::new(0.2, 0.1 * j as f64, -0.4 * k as f64))
                .collect(),
        })
        .collect();
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points).unwrap();
    let controls = adjoint.control_cotangents(&columns).unwrap();
    let nodes = adjoint.mesh_node_cotangents(&input, &columns).unwrap();
    let angles = adjoint.flight_condition_cotangents(&columns).unwrap();
    for k in 0..2 {
        for (name, derivative) in &controls[k] {
            let mut plus = condition.clone();
            let mut minus = condition.clone();
            plus.control_deflections
                .iter_mut()
                .find(|(n, _)| n == name)
                .unwrap()
                .1 += 1e-6;
            minus
                .control_deflections
                .iter_mut()
                .find(|(n, _)| n == name)
                .unwrap()
                .1 -= 1e-6;
            close(
                *derivative,
                (functional(&solve(&input, &plus).1, &points, &columns[k])
                    - functional(&solve(&input, &minus).1, &points, &columns[k]))
                    / 2e-6,
            );
        }
        for alpha in [true, false] {
            let mut plus = condition.clone();
            let mut minus = condition.clone();
            if alpha {
                plus.alpha += 1e-6;
                minus.alpha -= 1e-6;
            } else {
                plus.beta += 1e-6;
                minus.beta -= 1e-6;
            }
            close(
                if alpha {
                    angles[k].by_alpha
                } else {
                    angles[k].by_beta
                },
                (functional(&solve(&input, &plus).1, &points, &columns[k])
                    - functional(&solve(&input, &minus).1, &points, &columns[k]))
                    / 2e-6,
            );
        }
        let mut plus = input.clone();
        let mut minus = input.clone();
        let mut analytic = 0.0;
        for (j, seed) in nodes[k].surfaces[0].iter().enumerate() {
            let direction = Vec3::new((j as f64).sin(), 0.2, (0.3 * j as f64).cos());
            plus.surfaces[0].nodes[j] += 1e-6 * direction;
            minus.surfaces[0].nodes[j] -= 1e-6 * direction;
            analytic += seed.dot(&direction);
        }
        close(
            analytic,
            (functional(&solve(&plus, &condition).1, &points, &columns[k])
                - functional(&solve(&minus, &condition).1, &points, &columns[k]))
                / 2e-6,
        );
    }
    assert!(adjoint.control_cotangents(&[]).unwrap().is_empty());
    assert!(
        adjoint
            .control_cotangents(&[StripLoadCotangent::zeros(1)])
            .is_err()
    );
}

#[test]
fn invalid_controls_and_stale_derivatives_are_rejected() {
    let input = mesh(true, false);
    for mode in 0..6 {
        let mut invalid = input.clone();
        let control = &mut invalid.surfaces[0].controls[0];
        match mode {
            0 => control.name.clear(),
            1 => control.hinge.pop().map(|_| ()).unwrap(),
            2 => control.hinge[0] = f64::NAN,
            3 => control.span = 2..8,
            4 => control.gain = f64::INFINITY,
            _ => control.hinge[0] = 4.0,
        }
        assert!(VortexLatticeSolver::new().set_mesh(&invalid).is_err());
    }
    let condition = flight();
    let (mut solver, _) = solve(&input, &condition);
    let mut stale = condition.clone();
    stale.alpha += 0.1;
    assert!(
        solver
            .control_circulation_derivative(&stale, "elevator")
            .is_err()
    );
    assert!(
        solver
            .control_circulation_derivative(&condition, "typo")
            .is_err()
    );
    stale.control_deflections.push(("typo".into(), 0.0));
    assert!(solver.solve(&stale).is_err());
}

#[test]
fn section_geometry_gradients_include_deflected_hinge_axes_and_coverage() {
    use pegasus_vlm::geometry::{Aircraft, ControlSurfaceDefinition, Section, Surface};
    use pegasus_vlm::gradients::GeometryVariable;
    for inferred_axis in [false, true] {
        let mut wing = Surface::new("Wing");
        wing.set_chordwise_panel_count(3);
        wing.set_spanwise_panel_count(3);
        wing.set_y_duplicate(0.0);
        for (j, (chord, hinge)) in [(1.2, 0.56), (0.7, 0.69)].into_iter().enumerate() {
            let mut section = Section::new(
                Vec3::new(0.3 * j as f64, 0.2 + 3.0 * j as f64, 0.2 * j as f64),
                chord,
                2.0 - 4.0 * j as f64,
            );
            section.add_control_surface(ControlSurfaceDefinition {
                name: "elevator".into(),
                hinge_x_fraction: hinge,
                hinge_vector: if inferred_axis {
                    Vec3::zeros()
                } else {
                    Vec3::new(0.1, 1.0, 0.2)
                },
                sign: -1.0,
                gain: 0.8,
            });
            wing.add_section(section);
        }
        let mut aircraft = Aircraft::new("Controlled wing");
        aircraft.add_surface(wing);
        aircraft.set_reference_area(6.0);
        aircraft.set_reference_chord(1.0);
        aircraft.set_reference_span(6.0);
        let condition = FlightCondition {
            control_deflections: vec![("elevator".into(), 0.14)],
            ..flight()
        };
        let mut calculator = GradientCalculator::new();
        calculator.set_aircraft(&aircraft);
        let gradients = calculator.compute_geometry_gradients(&condition).unwrap();
        let evaluate = |a: &Aircraft| {
            let mut solver = VortexLatticeSolver::new();
            solver.set_aircraft(a).unwrap();
            let solution = solver.solve(&condition).unwrap();
            values(&ForceCalculator::new().compute_for_aircraft(
                a,
                solver.mesh(),
                &solution,
                &condition,
            ))
        };
        for (j, gradient) in gradients.sections.iter().enumerate() {
            for variable in [
                GeometryVariable::Chord,
                GeometryVariable::Twist,
                GeometryVariable::LeadingEdgeX,
                GeometryVariable::LeadingEdgeY,
                GeometryVariable::LeadingEdgeZ,
            ] {
                let perturb = |h: f64| {
                    let mut a = aircraft.clone();
                    let s = &mut a.surfaces_mut()[0].sections_mut()[j];
                    match variable {
                        GeometryVariable::Chord => s.set_chord(s.chord() + h),
                        GeometryVariable::Twist => s.set_twist_degrees(s.twist_degrees() + h),
                        _ => {
                            let mut p = *s.leading_edge();
                            let d = match variable {
                                GeometryVariable::LeadingEdgeX => 0,
                                GeometryVariable::LeadingEdgeY => 1,
                                _ => 2,
                            };
                            p[d] += h;
                            s.set_leading_edge(p);
                        }
                    }
                    evaluate(&a)
                };
                let p = perturb(1e-6);
                let m = perturb(-1e-6);
                for (k, output) in OUTPUTS.iter().enumerate() {
                    close(gradient.get(variable).get(*output), (p[k] - m[k]) / 2e-6);
                }
            }
        }
    }
}

#[test]
fn leading_edge_controls_and_shared_commands_on_multiple_surfaces() {
    let mut input = mesh(true, false);
    input.surfaces[0].controls[0].leading_edge = true;
    input.surfaces[0].controls[0].hinge = vec![0.6; 5];
    let mut second = input.surfaces[0].clone();
    second.name = "Tail".into();
    for point in &mut second.nodes {
        *point += Vec3::new(3.0, 0.0, 0.5);
    }
    second.controls[0].gain = -0.4;
    input.surfaces.push(second);
    let condition = flight();
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&input).unwrap();
    let gradient = calculator.compute_control_gradients(&condition).unwrap();
    assert_eq!(gradient.len(), 2);
    let mut plus = condition.clone();
    let mut minus = condition.clone();
    plus.control_deflections[0].1 += 1e-6;
    minus.control_deflections[0].1 -= 1e-6;
    let p = values(&solve(&input, &plus).1);
    let m = values(&solve(&input, &minus).1);
    let elevator = gradient.iter().find(|g| g.name == "elevator").unwrap();
    for (k, output) in OUTPUTS.iter().enumerate() {
        close(elevator.by_deflection.get(*output), (p[k] - m[k]) / 2e-6);
    }
    calculator.set_output_selection(OutputSelection::none());
    for g in calculator.compute_control_gradients(&condition).unwrap() {
        for output in OUTPUTS {
            assert_eq!(g.by_deflection.get(output), 0.0);
        }
    }
}

#[test]
fn cached_deflected_solve_matches_fresh_solves_after_condition_changes() {
    let input = mesh(true, true);
    let mut condition = flight();
    let (mut solver, _) = solve(&input, &condition);
    for mode in 0..5 {
        match mode {
            0 => condition.alpha += 0.03,
            1 => condition.pitch_rate -= 0.01,
            2 => condition.control_deflections[0].1 += 0.03,
            3 => condition.mach += 0.1,
            _ => condition.control_deflections.clear(),
        }
        let cached = solver.solve(&condition).unwrap();
        let (fresh, _) = solve(&input, &condition);
        assert!((&cached.circulation - &fresh.solution().circulation).norm() < 1e-12);
    }
}

#[test]
fn mirrored_controls_produce_symmetric_and_antisymmetric_responses() {
    let mut input = mesh(true, true);
    input.reference.moment_reference_point.y = 0.0;
    let surface = &mut input.surfaces[0];
    for point in &mut surface.nodes {
        point.z = 0.0;
    }
    surface.controls.truncate(1);
    let condition = FlightCondition {
        control_deflections: vec![("elevator".into(), 0.1)],
        ..FlightCondition::default()
    };
    let symmetric = solve(&input, &condition).1.coefficients;
    assert!(symmetric.lift.abs() > 0.01);
    assert!(symmetric.roll.abs() < 1e-12);
    input.surfaces[0].controls[0].duplicate_sign = -1.0;
    let antisymmetric = solve(&input, &condition).1.coefficients;
    assert!(antisymmetric.lift.abs() < 1e-12);
    assert!(antisymmetric.roll.abs() > 0.01);
}
