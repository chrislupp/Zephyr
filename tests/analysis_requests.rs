use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::analysis::{
    Analysis, AnalysisRequest, Details, GradientRequest, GradientValue, Input, Output, SolverMode,
};
use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{MeshInput, SurfaceMesh};
use pegasus_vlm::gradients::{AeroOutput, GradientCalculator};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn mesh(nx: usize, ny: usize) -> MeshInput {
    let mut s = SurfaceMesh::new("wing", nx, ny);
    s.y_duplicated = true;
    for i in 0..=nx {
        for j in 0..=ny {
            let eta = j as f64 / ny as f64;
            let chord = i as f64 / nx as f64;
            s.set_node(
                i,
                j,
                Vec3::new(
                    0.7 * eta + chord * (1.0 - 0.3 * eta),
                    4.0 * eta,
                    0.4 * eta - 0.08 * chord * eta,
                ),
            )
            .unwrap();
        }
    }
    let mut m = MeshInput::default();
    m.surfaces.push(s);
    m.reference.area = 6.8;
    m.reference.chord = 0.9;
    m.reference.span = 8.0;
    m.reference.moment_reference_point = Vec3::new(0.25, 0.3, 0.12);
    m
}
fn condition() -> FlightCondition {
    FlightCondition {
        alpha: 0.09,
        mach: 0.25,
        pitch_rate: 0.002,
        ..FlightCondition::default()
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-10 + 1e-8 * b.abs(),
        "{a:.15e} != {b:.15e}; error {}",
        (a - b).abs()
    );
}
fn scalar(v: &GradientValue) -> f64 {
    match v {
        GradientValue::Scalar(v) => *v,
        _ => panic!("expected scalar"),
    }
}
fn value(m: &MeshInput, c: &FlightCondition, mode: SolverMode, o: Output) -> f64 {
    Analysis::from_mesh(m)
        .evaluate(
            c,
            &AnalysisRequest {
                mode,
                outputs: vec![o],
                ..AnalysisRequest::default()
            },
        )
        .unwrap()
        .outputs[&o]
}

#[test]
fn selective_outputs_match_legacy_and_half_wing() {
    for (nx, ny) in [(1, 3), (3, 5), (4, 8)] {
        let m = mesh(nx, ny);
        let c = condition();
        let mut solver = VortexLatticeSolver::new();
        solver.set_mesh(&m).unwrap();
        let solution = solver.solve(&c).unwrap();
        let full = ForceCalculator::new().compute_with_references(
            &m.reference,
            solver.mesh(),
            &solution,
            &c,
        );
        let expected = [
            full.coefficients.lift,
            full.coefficients.drag,
            full.coefficients.side,
            full.coefficients.roll,
            full.coefficients.pitch,
            full.coefficients.yaw,
            full.near_field_induced_drag,
            full.viscous_drag,
            full.far_field_lift,
            full.far_field_side,
            full.far_field_induced_drag,
            full.span_efficiency,
        ];
        for mode in [SolverMode::Full, SolverMode::HalfWing] {
            let mut analysis = Analysis::from_mesh(&m);
            let request = AnalysisRequest {
                mode,
                outputs: Output::ALL.to_vec(),
                details: Details {
                    circulation: true,
                    induced_velocity: true,
                    surface_loads: true,
                },
                ..AnalysisRequest::default()
            };
            let result = analysis.evaluate(&c, &request).unwrap();
            for (o, e) in Output::ALL.iter().zip(expected) {
                close(result.outputs[o], e);
            }
            for (a, b) in result
                .circulation
                .as_ref()
                .unwrap()
                .iter()
                .zip(solution.circulation.iter())
            {
                close(*a, *b);
            }
            for (a, b) in result
                .induced_velocity
                .as_ref()
                .unwrap()
                .iter()
                .zip(solution.induced_velocity.iter())
            {
                close(*a, *b);
            }
            assert_eq!(
                result.surface_loads.as_ref().unwrap().len(),
                full.surface_loads.len()
            );
            for (a, b) in result
                .surface_loads
                .as_ref()
                .unwrap()
                .iter()
                .zip(&full.surface_loads)
            {
                assert_eq!(a.name, b.name);
                close(a.coefficients.lift, b.coefficients.lift);
                close(a.coefficients.pitch, b.coefficients.pitch);
                for (a, b) in a.strip_loads.iter().zip(&b.strip_loads) {
                    assert_eq!(a.strip_index, b.strip_index);
                    close(a.coefficients.lift, b.coefficients.lift);
                }
            }
            assert_eq!(
                result.work.unknowns,
                nx * ny * if mode == SolverMode::HalfWing { 1 } else { 2 }
            );
            assert_eq!(result.work.primal_solves, 1);
            assert_eq!(result.work.factorizations, 1);
        }
    }
}

#[test]
fn requests_skip_work_and_extend_caches() {
    let m = mesh(3, 6);
    let c = condition();
    let mut a = Analysis::from_mesh(&m);
    let empty = a.evaluate(&c, &AnalysisRequest::default()).unwrap();
    assert_eq!(empty.work.unknowns, 0);
    assert!(empty.outputs.is_empty());
    let mut r = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::FarFieldLift],
        ..AnalysisRequest::default()
    };
    let far = a.evaluate(&c, &r).unwrap();
    assert_eq!(far.work.near_field_panels, 0);
    assert_eq!(far.work.velocity_products, 0);
    assert_eq!(far.work.trefftz_interactions, 0);
    assert_eq!(far.work.influence_evaluations, 2 * 18 * 18);
    assert_eq!(far.work.operator_bytes, 8 * 4 * 18 * 18);
    assert!(far.circulation.is_none());
    assert!(far.surface_loads.is_none());
    r.outputs = vec![Output::Lift];
    let near = a.evaluate(&c, &r).unwrap();
    assert_eq!(near.outputs.len(), 1);
    assert_eq!(near.work.factorizations, 0);
    assert_eq!(near.work.primal_solves, 0);
    assert_eq!(near.work.velocity_products, 1);
    assert_eq!(near.work.operator_bytes, 8 * 7 * 18 * 18);
    assert_eq!(near.work.moment_panels, 0);
    assert_eq!(near.work.trefftz_interactions, 0);
    let hit = a.evaluate(&c, &r).unwrap();
    assert_eq!(hit.work.near_field_panels, 0);
    assert_eq!(hit.work.influence_evaluations, 0);
    assert_eq!(hit.outputs, near.outputs);
    r.outputs.clear();
    r.gradients = vec![GradientRequest {
        output: Output::Lift,
        input: Input::Alpha,
    }];
    let grad = a.evaluate(&c, &r).unwrap();
    assert!(grad.outputs.is_empty());
    assert_eq!(grad.gradients.len(), 1);
    assert_eq!(grad.work.adjoint_solves, 1);
    assert_eq!(grad.work.mach_interactions, 0);
    assert_eq!(grad.work.geometry_interactions, 0);
    assert_eq!(grad.work.factorizations, 0);
    assert_eq!(grad.work.primal_solves, 0);
    let hit = a.evaluate(&c, &r).unwrap();
    assert_eq!(hit.work.adjoint_solves, 0);
    r.gradients.push(GradientRequest {
        output: Output::Lift,
        input: Input::Mach,
    });
    let extended = a.evaluate(&c, &r).unwrap();
    assert_eq!(extended.work.adjoint_solves, 0);
    assert!(extended.work.mach_interactions > 0);
    assert_eq!(extended.work.alpha_derivative_panels, 0);
}

#[test]
fn selected_flight_gradients_match_full_and_finite_differences() {
    let m = mesh(3, 5);
    let c = condition();
    let inputs = [
        Input::Alpha,
        Input::Mach,
        Input::PitchRate,
        Input::MomentReference(0),
        Input::MomentReference(1),
        Input::MomentReference(2),
    ];
    for o in Output::ALL {
        let request = AnalysisRequest {
            gradients: inputs
                .iter()
                .cloned()
                .map(|input| GradientRequest { output: o, input })
                .collect(),
            ..AnalysisRequest::default()
        };
        let full = Analysis::from_mesh(&m).evaluate(&c, &request).unwrap();
        let half = Analysis::from_mesh(&m)
            .evaluate(
                &c,
                &AnalysisRequest {
                    mode: SolverMode::HalfWing,
                    ..request
                },
            )
            .unwrap();
        for (i, input) in inputs.iter().enumerate() {
            let a = scalar(&full.gradients[i].value);
            let b = scalar(&half.gradients[i].value);
            close(a, b);
            let h = 1e-5;
            let mut plus = c.clone();
            let mut minus = c.clone();
            let mut pm = m.clone();
            let mut mm = m.clone();
            match input {
                Input::Alpha => {
                    plus.alpha += h;
                    minus.alpha -= h;
                }
                Input::Mach => {
                    plus.mach += h;
                    minus.mach -= h;
                }
                Input::PitchRate => {
                    plus.pitch_rate += h;
                    minus.pitch_rate -= h;
                }
                Input::MomentReference(k) => {
                    pm.reference.moment_reference_point[*k] += h;
                    mm.reference.moment_reference_point[*k] -= h;
                }
                _ => unreachable!(),
            }
            let fd = (value(&pm, &plus, SolverMode::HalfWing, o)
                - value(&mm, &minus, SolverMode::HalfWing, o))
                / (2.0 * h);
            assert!(
                (fd - b).abs() < 2e-7 + 2e-6 * fd.abs(),
                "{o:?}/{input:?}: analytic {b}, finite difference {fd}"
            );
        }
    }
}

#[test]
fn geometry_gradients_fold_mirrored_inputs() {
    let m = mesh(2, 4);
    let c = condition();
    for o in [
        Output::Lift,
        Output::Drag,
        Output::Roll,
        Output::Pitch,
        Output::Yaw,
        Output::FarFieldLift,
        Output::FarFieldInducedDrag,
        Output::SpanEfficiency,
    ] {
        let req = AnalysisRequest {
            gradients: vec![GradientRequest {
                output: o,
                input: Input::MeshNodes(0),
            }],
            ..AnalysisRequest::default()
        };
        let full = Analysis::from_mesh(&m).evaluate(&c, &req).unwrap();
        let half = Analysis::from_mesh(&m)
            .evaluate(
                &c,
                &AnalysisRequest {
                    mode: SolverMode::HalfWing,
                    ..req
                },
            )
            .unwrap();
        let (GradientValue::MeshNodes(f), GradientValue::MeshNodes(h)) =
            (&full.gradients[0].value, &half.gradients[0].value)
        else {
            panic!()
        };
        for (a, b) in f.iter().zip(h) {
            for k in 0..3 {
                close(a[k], b[k]);
            }
        }
        // Perturb interior and tip nodes; y=0 root positions stay on the plane.
        let mut pm = m.clone();
        let mut mm = m.clone();
        let step = 1e-7;
        let mut analytic = 0.0;
        for (i, node) in m.surfaces[0].nodes.iter().enumerate() {
            let direction = Vec3::new(
                (i as f64 * 0.7).sin(),
                if node.y == 0.0 {
                    0.0
                } else {
                    (i as f64 * 0.4).cos()
                },
                (i as f64 * 0.9).cos(),
            );
            pm.surfaces[0].nodes[i] += step * direction;
            mm.surfaces[0].nodes[i] -= step * direction;
            analytic += h[i].dot(&direction);
        }
        let fd = (value(&pm, &c, SolverMode::HalfWing, o)
            - value(&mm, &c, SolverMode::HalfWing, o))
            / (2.0 * step);
        assert!(
            (fd - analytic).abs() < 2e-7 + 2e-6 * fd.abs(),
            "{o:?}: analytic {analytic}, fd {fd}"
        );
    }
}

#[test]
fn gradients_match_existing_calculator() {
    let m = mesh(3, 5);
    let c = condition();
    let mut old = GradientCalculator::new();
    old.set_mesh(&m).unwrap();
    let expected = old.compute_flight_condition_gradients(&c).unwrap();
    for (o, ao) in [
        (Output::Lift, AeroOutput::Lift),
        (Output::Drag, AeroOutput::Drag),
        (Output::Pitch, AeroOutput::Pitch),
        (Output::FarFieldInducedDrag, AeroOutput::InducedDrag),
    ] {
        for (input, e) in [
            (Input::Alpha, expected.by_alpha.get(ao)),
            (Input::Beta, expected.by_beta.get(ao)),
            (Input::Mach, expected.by_mach.get(ao)),
            (Input::PitchRate, expected.by_pitch_rate.get(ao)),
        ] {
            let result = Analysis::from_mesh(&m)
                .evaluate(
                    &c,
                    &AnalysisRequest {
                        gradients: vec![GradientRequest { output: o, input }],
                        ..AnalysisRequest::default()
                    },
                )
                .unwrap();
            close(scalar(&result.gradients[0].value), e);
        }
    }
}

#[test]
fn restrictions_do_not_affect_full_mode_and_failed_requests_do_not_poison_cache() {
    let m = mesh(2, 4);
    let mut a = Analysis::from_mesh(&m);
    let c = condition();
    let mut req = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift],
        ..AnalysisRequest::default()
    };
    let expected = a.evaluate(&c, &req).unwrap().outputs[&Output::Lift];
    let lateral = FlightCondition {
        beta: 0.05,
        roll_rate: 0.01,
        ..c.clone()
    };
    assert!(a.evaluate(&lateral, &req).is_err());
    close(
        a.evaluate(&c, &req).unwrap().outputs[&Output::Lift],
        expected,
    );
    req.gradients.push(GradientRequest {
        output: Output::Roll,
        input: Input::RollRate,
    });
    assert!(a.evaluate(&c, &req).is_err());
    req.mode = SolverMode::Full;
    assert!(a.evaluate(&lateral, &req).is_ok());
    req.gradients.clear();
    req.mode = SolverMode::HalfWing;
    close(
        a.evaluate(&c, &req).unwrap().outputs[&Output::Lift],
        expected,
    );
    let mut changed = m.clone();
    changed.surfaces[0].nodes[4].x += 0.1;
    a.set_mesh(&changed);
    close(
        a.evaluate(&c, &req).unwrap().outputs[&Output::Lift],
        value(&changed, &c, SolverMode::HalfWing, Output::Lift),
    );
}

mod common;

#[test]
fn aircraft_sections_controls_and_viscous_loads_preserve_behavior() {
    use pegasus_vlm::geometry::ControlSurfaceDefinition;
    use pegasus_vlm::gradients::GeometryVariable;
    let mut aircraft = common::tapered_twisted_wing_mirrored();
    aircraft.set_moment_reference_point(Vec3::new(0.3, 0.2, 0.1));
    let c = condition();
    for variable in [
        GeometryVariable::Chord,
        GeometryVariable::Twist,
        GeometryVariable::LeadingEdgeX,
        GeometryVariable::LeadingEdgeY,
        GeometryVariable::LeadingEdgeZ,
    ] {
        let req = AnalysisRequest {
            outputs: Output::ALL.to_vec(),
            gradients: [Output::Lift, Output::Pitch, Output::FarFieldInducedDrag]
                .into_iter()
                .map(|output| GradientRequest {
                    output,
                    input: Input::Section {
                        surface: 0,
                        section: 1,
                        variable,
                    },
                })
                .collect(),
            ..AnalysisRequest::default()
        };
        let full = Analysis::from_aircraft(&aircraft)
            .evaluate(&c, &req)
            .unwrap();
        let half = Analysis::from_aircraft(&aircraft)
            .evaluate(
                &c,
                &AnalysisRequest {
                    mode: SolverMode::HalfWing,
                    ..req
                },
            )
            .unwrap();
        for o in Output::ALL {
            close(full.outputs[&o], half.outputs[&o]);
        }
        for (a, b) in full.gradients.iter().zip(&half.gradients) {
            close(scalar(&a.value), scalar(&b.value));
        }
    }
    for s in aircraft.surfaces_mut()[0].sections_mut() {
        s.add_control_surface(ControlSurfaceDefinition {
            name: "flap".into(),
            hinge_x_fraction: 0.7,
            hinge_vector: Vec3::y(),
            sign: 1.0,
            gain: 1.0,
        });
        s.set_viscous_polar([-1.0, 0.03, 0.0, 0.01, 1.0, 0.03]);
    }
    aircraft.set_profile_drag_coefficient(0.015);
    let controlled = FlightCondition {
        control_deflections: vec![("flap".into(), 0.07)],
        ..c.clone()
    };
    for mode in [SolverMode::Full, SolverMode::HalfWing] {
        let req = AnalysisRequest {
            mode,
            outputs: Output::ALL.to_vec(),
            gradients: [
                Output::Drag,
                Output::Pitch,
                Output::Roll,
                Output::ViscousDrag,
            ]
            .into_iter()
            .flat_map(|output| {
                [Input::Alpha, Input::Mach, Input::Control("flap".into())]
                    .into_iter()
                    .map(move |input| GradientRequest { output, input })
            })
            .collect(),
            ..AnalysisRequest::default()
        };
        let result = Analysis::from_aircraft(&aircraft)
            .evaluate(&controlled, &req)
            .unwrap();
        let mut solver = VortexLatticeSolver::new();
        solver.set_aircraft(&aircraft).unwrap();
        let solution = solver.solve(&controlled).unwrap();
        let old = ForceCalculator::new().compute_for_aircraft(
            &aircraft,
            solver.mesh(),
            &solution,
            &controlled,
        );
        close(result.outputs[&Output::Drag], old.coefficients.drag);
        close(result.outputs[&Output::Pitch], old.coefficients.pitch);
        close(result.outputs[&Output::ViscousDrag], old.viscous_drag);
        for g in &result.gradients {
            let mut plus = controlled.clone();
            let mut minus = controlled.clone();
            let h = 1e-5;
            match g.request.input {
                Input::Alpha => {
                    plus.alpha += h;
                    minus.alpha -= h;
                }
                Input::Mach => {
                    plus.mach += h;
                    minus.mach -= h;
                }
                Input::Control(_) => {
                    plus.control_deflections[0].1 += h;
                    minus.control_deflections[0].1 -= h;
                }
                _ => unreachable!(),
            }
            let r = AnalysisRequest {
                mode,
                outputs: vec![g.request.output],
                ..AnalysisRequest::default()
            };
            let fd = (Analysis::from_aircraft(&aircraft)
                .evaluate(&plus, &r)
                .unwrap()
                .outputs[&g.request.output]
                - Analysis::from_aircraft(&aircraft)
                    .evaluate(&minus, &r)
                    .unwrap()
                    .outputs[&g.request.output])
                / (2.0 * h);
            assert!(
                (fd - scalar(&g.value)).abs() < 2e-7 + 2e-6 * fd.abs(),
                "{:?}: {} != {fd}",
                g.request,
                scalar(&g.value)
            );
        }
    }
}

#[test]
fn multiple_components_and_left_half_meshes_match() {
    let mut m = mesh(2, 5);
    let mut tail = m.surfaces[0].clone();
    tail.name = "tail".into();
    tail.component_index = 1;
    for node in &mut tail.nodes {
        node.x = 3.0 + 0.5 * node.x;
        node.y *= 0.7;
        node.z += 0.5;
    }
    m.surfaces.push(tail);
    for left in [false, true] {
        let mut input = m.clone();
        if left {
            for s in &mut input.surfaces {
                for p in &mut s.nodes {
                    p.y = -p.y;
                }
            }
        }
        let req = AnalysisRequest {
            outputs: Output::ALL.to_vec(),
            gradients: vec![
                GradientRequest {
                    output: Output::FarFieldInducedDrag,
                    input: Input::MeshNodes(1),
                },
                GradientRequest {
                    output: Output::Lift,
                    input: Input::MeshNodes(0),
                },
            ],
            ..AnalysisRequest::default()
        };
        let full = Analysis::from_mesh(&input)
            .evaluate(&condition(), &req)
            .unwrap();
        let half = Analysis::from_mesh(&input)
            .evaluate(
                &condition(),
                &AnalysisRequest {
                    mode: SolverMode::HalfWing,
                    ..req
                },
            )
            .unwrap();
        for o in Output::ALL {
            close(full.outputs[&o], half.outputs[&o]);
        }
        for (a, b) in full.gradients.iter().zip(&half.gradients) {
            let (GradientValue::MeshNodes(a), GradientValue::MeshNodes(b)) = (&a.value, &b.value)
            else {
                panic!()
            };
            for (a, b) in a.iter().zip(b) {
                for k in 0..3 {
                    close(a[k], b[k]);
                }
            }
        }
    }
}

#[test]
fn exact_zero_and_normalizing_derivatives_skip_adjoint_work() {
    let m = mesh(2, 4);
    let c = condition();
    let mut a = Analysis::from_mesh(&m);
    let zero = AnalysisRequest {
        gradients: vec![
            GradientRequest {
                output: Output::Lift,
                input: Input::DynamicPressure,
            },
            GradientRequest {
                output: Output::Lift,
                input: Input::ReferenceChord,
            },
            GradientRequest {
                output: Output::SpanEfficiency,
                input: Input::ReferenceArea,
            },
        ],
        ..AnalysisRequest::default()
    };
    let result = a.evaluate(&c, &zero).unwrap();
    assert_eq!(result.work.unknowns, 0);
    assert_eq!(result.gradients.len(), 3);
    for g in result.gradients {
        assert_eq!(scalar(&g.value), 0.0);
    }
    let r = AnalysisRequest {
        gradients: vec![GradientRequest {
            output: Output::Lift,
            input: Input::ReferenceArea,
        }],
        ..AnalysisRequest::default()
    };
    let result = a.evaluate(&c, &r).unwrap();
    assert_eq!(result.work.adjoint_solves, 0);
    close(
        scalar(&result.gradients[0].value),
        -value(&m, &c, SolverMode::Full, Output::Lift) / m.reference.area,
    );
    let r = AnalysisRequest {
        gradients: vec![GradientRequest {
            output: Output::Lift,
            input: Input::Mach,
        }],
        ..AnalysisRequest::default()
    };
    let result = a.evaluate(&c, &r).unwrap();
    assert_eq!(result.work.alpha_derivative_panels, 0);
    assert_eq!(result.work.beta_derivative_panels, 0);
    assert_eq!(result.work.adjoint_solves, 1);
    assert!(result.work.mach_interactions > 0);
}

#[test]
fn cache_updates_for_mach_rates_and_symmetric_controls() {
    use pegasus_vlm::geometry::ControlSurfaceDefinition;
    let mut aircraft = common::tapered_twisted_wing_mirrored();
    for section in aircraft.surfaces_mut()[0].sections_mut() {
        section.add_control_surface(ControlSurfaceDefinition {
            name: "flap".into(),
            hinge_x_fraction: 0.7,
            hinge_vector: Vec3::y(),
            sign: 1.0,
            gain: 1.0,
        });
    }
    let mut a = Analysis::from_aircraft(&aircraft);
    let req = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift],
        ..AnalysisRequest::default()
    };
    let mut c = condition();
    a.evaluate(&c, &req).unwrap();
    c.alpha += 0.02;
    let alpha = a.evaluate(&c, &req).unwrap();
    assert_eq!(alpha.work.factorizations, 0);
    assert_eq!(alpha.work.primal_solves, 1);
    c.pitch_rate += 0.001;
    let rate = a.evaluate(&c, &req).unwrap();
    assert_eq!(rate.work.factorizations, 0);
    assert_eq!(rate.work.primal_solves, 1);
    c.control_deflections = vec![("flap".into(), 0.07)];
    let control = a.evaluate(&c, &req).unwrap();
    assert_eq!(control.work.factorizations, 1);
    assert_eq!(control.work.influence_evaluations, 0);
    close(
        control.outputs[&Output::Lift],
        Analysis::from_aircraft(&aircraft)
            .evaluate(&c, &req)
            .unwrap()
            .outputs[&Output::Lift],
    );
    c.mach += 0.1;
    let mach = a.evaluate(&c, &req).unwrap();
    assert_eq!(mach.work.factorizations, 1);
    assert!(mach.work.influence_evaluations > 0);
    close(
        mach.outputs[&Output::Lift],
        Analysis::from_aircraft(&aircraft)
            .evaluate(&c, &req)
            .unwrap()
            .outputs[&Output::Lift],
    );
}

#[test]
fn invalid_half_geometries_are_rejected_before_assembly() {
    let c = condition();
    let req = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift],
        ..AnalysisRequest::default()
    };
    let m = mesh(2, 4);
    let mut cases = Vec::new();
    let mut unmatched = m.clone();
    unmatched.surfaces[0].y_duplicated = false;
    cases.push(unmatched);
    let mut offset = m.clone();
    offset.surfaces[0].y_duplicate_plane = 0.3;
    cases.push(offset);
    let mut crossing = m.clone();
    for p in &mut crossing.surfaces[0].nodes {
        p.y -= 2.0;
    }
    cases.push(crossing);
    let mut center = m.clone();
    for p in &mut center.surfaces[0].nodes {
        p.y = 0.0;
    }
    cases.push(center);
    for input in cases {
        assert!(Analysis::from_mesh(&input).evaluate(&c, &req).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, -0.1, 1.0] {
        let c = FlightCondition {
            mach: bad,
            ..condition()
        };
        assert!(Analysis::from_mesh(&m).evaluate(&c, &req).is_err());
    }
}

#[cfg(feature = "parallel")]
#[test]
fn selective_parallel_results_are_independent_of_thread_count() {
    use pegasus_vlm::analysis::AnalysisResult;
    fn bits(r: &AnalysisResult) -> Vec<u64> {
        let mut result: Vec<_> = r.outputs.values().map(|v| v.to_bits()).collect();
        for g in &r.gradients {
            match &g.value {
                GradientValue::Scalar(v) => result.push(v.to_bits()),
                GradientValue::MeshNodes(nodes) => {
                    result.extend(nodes.iter().flat_map(|v| v.iter().map(|x| x.to_bits())))
                }
            }
        }
        result
    }
    for (nx, ny) in [(2, 8), (6, 24)] {
        let m = mesh(nx, ny);
        for mode in [SolverMode::Full, SolverMode::HalfWing] {
            let request = AnalysisRequest {
                mode,
                outputs: vec![Output::Lift, Output::Drag, Output::Pitch],
                gradients: [Output::Lift, Output::Drag, Output::Pitch]
                    .into_iter()
                    .flat_map(|output| {
                        [Input::Alpha, Input::MeshNodes(0)]
                            .into_iter()
                            .map(move |input| GradientRequest { output, input })
                    })
                    .collect(),
                ..AnalysisRequest::default()
            };
            let mut reference = None;
            for threads in [1, 3, 10] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let actual = pool.install(|| {
                    let direct = Analysis::from_mesh(&m)
                        .evaluate(&condition(), &request)
                        .unwrap();
                    // Exercise AIC-only then velocity-only assembly as well.
                    let mut staged = Analysis::from_mesh(&m);
                    let far = staged
                        .evaluate(
                            &condition(),
                            &AnalysisRequest {
                                mode,
                                outputs: vec![Output::FarFieldInducedDrag],
                                ..AnalysisRequest::default()
                            },
                        )
                        .unwrap();
                    let extended = staged.evaluate(&condition(), &request).unwrap();
                    assert_eq!(bits(&direct), bits(&extended));
                    assert_eq!(
                        direct.work.unknowns,
                        nx * ny * if mode == SolverMode::HalfWing { 1 } else { 2 }
                    );
                    let mut values = bits(&far);
                    values.extend(bits(&direct));
                    values
                });
                if let Some(expected) = &reference {
                    assert_eq!(&actual, expected);
                } else {
                    reference = Some(actual);
                }
            }
        }
    }
}

#[test]
fn half_wing_mesh_gradient_above_parallel_threshold_matches_difference() {
    let mut m = mesh(6, 24);
    let c = condition();
    let output = Output::Drag;
    let request = AnalysisRequest {
        mode: SolverMode::HalfWing,
        gradients: vec![GradientRequest {
            output,
            input: Input::MeshNodes(0),
        }],
        ..AnalysisRequest::default()
    };
    let result = Analysis::from_mesh(&m).evaluate(&c, &request).unwrap();
    let GradientValue::MeshNodes(nodes) = &result.gradients[0].value else {
        panic!()
    };
    let index = 3 * 25 + 12;
    let step = 1e-7;
    m.surfaces[0].nodes[index].z += step;
    let plus = value(&m, &c, SolverMode::HalfWing, output);
    m.surfaces[0].nodes[index].z -= 2.0 * step;
    let minus = value(&m, &c, SolverMode::HalfWing, output);
    close(nodes[index].z, (plus - minus) / (2.0 * step));
}
