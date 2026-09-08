mod common;

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3, VecX};
use pegasus_vlm::gradients::{
    AeroOutput, GeometryVariable, GradientCalculator, OutputSelection, TotalDerivative,
    aic_element_gradient, horseshoe_velocity_pg_derivative, horseshoe_velocity_vjp_bound_left,
    horseshoe_velocity_vjp_bound_right, horseshoe_velocity_vjp_field_point,
    trefftz_circulation_gradients,
};
use pegasus_vlm::solver::{
    FlightCondition, SolutionState, VortexLatticeSolver, horseshoe_vortex_velocity,
};

use common::{rectangular_mesh, rectangular_wing};

fn coefficient(
    aircraft: &pegasus_vlm::geometry::Aircraft,
    condition: &FlightCondition,
    output: AeroOutput,
) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(aircraft).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result =
        ForceCalculator::new().compute_for_aircraft(aircraft, solver.mesh(), &solution, condition);
    match output {
        AeroOutput::Lift => result.coefficients.lift,
        AeroOutput::Drag => result.coefficients.drag,
        AeroOutput::Side => result.coefficients.side,
        AeroOutput::Roll => result.coefficients.roll,
        AeroOutput::Pitch => result.coefficients.pitch,
        AeroOutput::Yaw => result.coefficients.yaw,
        AeroOutput::InducedDrag => result.far_field_induced_drag,
    }
}

fn mesh_coefficient(
    mesh: &pegasus_vlm::geometry::MeshInput,
    condition: &FlightCondition,
    output: AeroOutput,
) -> f64 {
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(mesh).unwrap();
    let solution = solver.solve(condition).unwrap();
    let result = ForceCalculator::new().compute_with_references(
        &mesh.reference,
        solver.mesh(),
        &solution,
        condition,
    );
    match output {
        AeroOutput::Lift => result.coefficients.lift,
        AeroOutput::Drag => result.coefficients.drag,
        AeroOutput::Side => result.coefficients.side,
        AeroOutput::Roll => result.coefficients.roll,
        AeroOutput::Pitch => result.coefficients.pitch,
        AeroOutput::Yaw => result.coefficients.yaw,
        AeroOutput::InducedDrag => result.far_field_induced_drag,
    }
}

#[test]
fn biot_savart_vjps_match_centered_differences() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let seed = Vec3::new(0.4, -0.7, 1.2);
    let step = 1.0e-7;
    let analytical = [
        horseshoe_velocity_vjp_field_point(&field, &left, &right, 0.92, 0.01, true, &seed),
        horseshoe_velocity_vjp_bound_left(&field, &left, &right, 0.92, 0.01, true, &seed),
        horseshoe_velocity_vjp_bound_right(&field, &left, &right, 0.92, 0.01, true, &seed),
    ];
    for (which, base) in [field, left, right].into_iter().enumerate() {
        for component in 0..3 {
            let mut plus = [field, left, right];
            let mut minus = plus;
            plus[which] = base;
            minus[which] = base;
            plus[which][component] += step;
            minus[which][component] -= step;
            let value = |points: [Vec3; 3]| {
                horseshoe_vortex_velocity(&points[0], &points[1], &points[2], 0.92, 0.01, true)
                    .dot(&seed)
            };
            let finite = (value(plus) - value(minus)) / (2.0 * step);
            assert_abs_diff_eq!(analytical[which][component], finite, epsilon = 2.0e-7);
        }
    }
}

/// The same check with every length scaled down by 1e-4.
///
/// `horseshoe_vortex_velocity` tests collinearity with the dimensionless
/// `cross_sq > EPSILON * a_sq * b_sq`, but the VJP kernel used to compare the
/// length⁴ quantity `cross_sq` against an absolute `EPSILON`. Below roughly
/// millimetre scale that is false for ordinary panels, so the derivative
/// dropped the bound-leg term the primal still evaluated — silently, and with
/// no error reported. At this scale `cross_sq ≈ 1e-16` against an `EPSILON` of
/// `1e-12`, while the primal's threshold is `1e-28`, so the two disagreed by
/// the whole bound contribution: the bound_right x-derivative came out with the
/// wrong sign entirely.
///
/// The tolerance is relative because the velocity itself grows as the geometry
/// shrinks; at unit scale these derivatives are O(1), here they are O(1e6).
#[test]
fn biot_savart_vjps_match_centered_differences_at_small_scale() {
    let scale = 1.0e-4;
    let field = Vec3::new(0.72, 0.31, 0.18) * scale;
    let left = Vec3::new(0.25, -0.4, 0.02) * scale;
    let right = Vec3::new(0.28, 0.55, -0.01) * scale;
    let seed = Vec3::new(0.4, -0.7, 1.2);
    // The core radius is derived from panel size, so it scales with the mesh.
    let core = 0.01 * scale;
    let step = 1.0e-6 * scale;
    let analytical = [
        horseshoe_velocity_vjp_field_point(&field, &left, &right, 0.92, core, true, &seed),
        horseshoe_velocity_vjp_bound_left(&field, &left, &right, 0.92, core, true, &seed),
        horseshoe_velocity_vjp_bound_right(&field, &left, &right, 0.92, core, true, &seed),
    ];
    for (which, base) in [field, left, right].into_iter().enumerate() {
        for component in 0..3 {
            let mut plus = [field, left, right];
            let mut minus = plus;
            plus[which] = base;
            minus[which] = base;
            plus[which][component] += step;
            minus[which][component] -= step;
            let value = |points: [Vec3; 3]| {
                horseshoe_vortex_velocity(&points[0], &points[1], &points[2], 0.92, core, true)
                    .dot(&seed)
            };
            let finite = (value(plus) - value(minus)) / (2.0 * step);
            let analytical = analytical[which][component];
            let magnitude = analytical.abs().max(finite.abs());
            assert!(
                (analytical - finite).abs() <= 1.0e-5 * magnitude,
                "vjp {which} component {component}: analytical {analytical:e} \
                 against finite difference {finite:e}"
            );
        }
    }
}

#[test]
fn biot_savart_pg_derivative_matches_centered_difference() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let pg = 0.91;
    let step = 1.0e-7;
    let analytical = horseshoe_velocity_pg_derivative(&field, &left, &right, pg, 0.01, true);
    let finite = (horseshoe_vortex_velocity(&field, &left, &right, pg + step, 0.01, true)
        - horseshoe_vortex_velocity(&field, &left, &right, pg - step, 0.01, true))
        / (2.0 * step);
    assert_abs_diff_eq!(analytical, finite, epsilon = 3.0e-8);
}

#[test]
fn aic_element_gradient_matches_centered_differences() {
    let field = Vec3::new(0.72, 0.31, 0.18);
    let left = Vec3::new(0.25, -0.4, 0.02);
    let right = Vec3::new(0.28, 0.55, -0.01);
    let normal = Vec3::new(0.1, -0.2, 0.97).normalize();
    let gradient = aic_element_gradient(&left, &right, &field, &normal, 0.02, 0.9, true);
    let step = 1.0e-7;
    for (which, analytical) in [
        gradient.bound_left,
        gradient.bound_right,
        gradient.control_point,
        gradient.normal,
    ]
    .into_iter()
    .enumerate()
    {
        for component in 0..3 {
            let mut plus = [left, right, field, normal];
            let mut minus = plus;
            plus[which][component] += step;
            minus[which][component] -= step;
            let value = |points: [Vec3; 4]| {
                horseshoe_vortex_velocity(&points[2], &points[0], &points[1], 0.9, 0.02, true)
                    .dot(&points[3])
            };
            let finite = (value(plus) - value(minus)) / (2.0 * step);
            assert_abs_diff_eq!(analytical[component], finite, epsilon = 3.0e-7);
        }
    }
}

#[test]
fn trefftz_circulation_gradient_matches_direct_perturbations() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    let solution = solver.solve(&condition).unwrap();
    let gradient =
        trefftz_circulation_gradients(solver.mesh(), &solution, 4.0, 0.01, 0.1).induced_drag;
    let step = 1.0e-7;
    for index in [
        0,
        solution.circulation.len() / 2,
        solution.circulation.len() - 1,
    ] {
        let mut plus = solution.clone();
        let mut minus = solution.clone();
        plus.circulation[index] += step;
        minus.circulation[index] -= step;
        let calculator = ForceCalculator::new();
        let plus_value = calculator.trefftz_plane(solver.mesh(), &plus, 4.0, 4.0).2;
        let minus_value = calculator.trefftz_plane(solver.mesh(), &minus, 4.0, 4.0).2;
        assert_abs_diff_eq!(
            gradient[index],
            (plus_value - minus_value) / (2.0 * step),
            epsilon = 1.0e-8
        );
    }
}

#[test]
fn adjoint_induced_drag_matches_full_resolves() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    let derivative =
        TotalDerivative::compute(&solver, &condition, AeroOutput::InducedDrag).unwrap();
    let step = 1.0e-6;
    let mut plus = condition.clone();
    let mut minus = condition.clone();
    plus.alpha += step;
    minus.alpha -= step;
    let finite_alpha = (coefficient(&aircraft, &plus, AeroOutput::InducedDrag)
        - coefficient(&aircraft, &minus, AeroOutput::InducedDrag))
        / (2.0 * step);
    assert_abs_diff_eq!(derivative.by_alpha, finite_alpha, epsilon = 2.0e-6);
}

#[test]
fn gradient_calculator_flight_derivatives_match_resolves() {
    let aircraft = rectangular_wing(3, 6);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 2.0 * DEGREES_TO_RADIANS,
        mach: 0.2,
        ..FlightCondition::default()
    };
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_flight_condition_gradients(&condition)
        .unwrap();
    let step = 1.0e-6;
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
        AeroOutput::InducedDrag,
    ] {
        let mut plus = condition.clone();
        let mut minus = condition.clone();
        plus.alpha += step;
        minus.alpha -= step;
        let finite = (coefficient(&aircraft, &plus, output)
            - coefficient(&aircraft, &minus, output))
            / (2.0 * step);
        assert_abs_diff_eq!(gradients.by_alpha.get(output), finite, epsilon = 3.0e-5);

        plus = condition.clone();
        minus = condition.clone();
        plus.mach += step;
        minus.mach -= step;
        let finite_mach = (coefficient(&aircraft, &plus, output)
            - coefficient(&aircraft, &minus, output))
            / (2.0 * step);
        assert_abs_diff_eq!(gradients.by_mach.get(output), finite_mach, epsilon = 3.0e-5);
    }
}

#[test]
fn output_selection_leaves_disabled_outputs_exactly_zero() {
    let aircraft = rectangular_wing(2, 4);
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Pitch));
    let gradients = calculator
        .compute_flight_condition_gradients(&FlightCondition::default())
        .unwrap();
    assert_eq!(gradients.by_alpha.coefficients.lift, 0.0);
    assert_eq!(gradients.by_alpha.coefficients.drag, 0.0);
    assert_eq!(gradients.by_alpha.induced_drag, 0.0);
    assert!(gradients.by_alpha.coefficients.pitch.is_finite());
}

#[test]
fn section_geometry_gradient_has_expected_shape_and_nonzero_chord_effect() {
    let aircraft = rectangular_wing(2, 4);
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Lift));
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let gradients = calculator.compute_geometry_gradients(&condition).unwrap();
    assert_eq!(gradients.sections.len(), 2);
    assert_eq!(gradients.design_variable_count(), 10);
    assert!(
        gradients.sections[0]
            .get(GeometryVariable::Chord)
            .coefficients
            .lift
            .abs()
            > 1.0e-6
    );
    assert_eq!(gradients.sections[0].by_chord.coefficients.drag, 0.0);
}

#[test]
fn direct_mesh_node_gradient_matches_external_resolve() {
    let mesh = rectangular_mesh(1, 2);
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&mesh).unwrap();
    calculator.set_output_selection(OutputSelection::only(AeroOutput::Lift));
    let condition = FlightCondition {
        alpha: 5.0 * DEGREES_TO_RADIANS,
        ..FlightCondition::default()
    };
    let gradients = calculator.compute_mesh_node_gradients(&condition).unwrap();
    assert_eq!(gradients.len(), 1);
    assert_eq!(gradients[0].lift.nodes.len(), 6);
    assert_eq!(gradients[0].drag.nodes, vec![Vec3::zeros(); 6]);

    let step = 1.0e-6;
    let mut plus = mesh.clone();
    let mut minus = mesh.clone();
    plus.surfaces[0].nodes[0].z += step;
    minus.surfaces[0].nodes[0].z -= step;
    let solve = |input: &pegasus_vlm::geometry::MeshInput| {
        let mut solver = VortexLatticeSolver::new();
        solver.set_mesh(input).unwrap();
        let solution = solver.solve(&condition).unwrap();
        ForceCalculator::new()
            .compute_with_references(&input.reference, solver.mesh(), &solution, &condition)
            .coefficients
            .lift
    };
    assert_abs_diff_eq!(
        gradients[0].lift.nodes[0].z,
        (solve(&plus) - solve(&minus)) / (2.0 * step),
        epsilon = 1.0e-7
    );
}

#[test]
fn direct_mesh_node_adjoint_matches_all_outputs() {
    let mesh = rectangular_mesh(2, 3);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        roll_rate: 0.03,
        pitch_rate: -0.02,
        yaw_rate: 0.01,
        ..FlightCondition::default()
    };
    let mut calculator = GradientCalculator::new();
    calculator.set_mesh(&mesh).unwrap();
    let gradients = calculator.compute_mesh_node_gradients(&condition).unwrap();
    let node_i = 1;
    let node_j = 1;
    let node_index = node_i * (mesh.surfaces[0].ny + 1) + node_j;
    let component = 2;
    let step = 1.0e-6;
    let mut plus = mesh.clone();
    let mut minus = mesh.clone();
    plus.surfaces[0].nodes[node_index][component] += step;
    minus.surfaces[0].nodes[node_index][component] -= step;
    for (output, analytical) in [
        (
            AeroOutput::Lift,
            gradients[0].lift.nodes[node_index][component],
        ),
        (
            AeroOutput::Drag,
            gradients[0].drag.nodes[node_index][component],
        ),
        (
            AeroOutput::Side,
            gradients[0].side.nodes[node_index][component],
        ),
        (
            AeroOutput::Roll,
            gradients[0].roll.nodes[node_index][component],
        ),
        (
            AeroOutput::Pitch,
            gradients[0].pitch.nodes[node_index][component],
        ),
        (
            AeroOutput::Yaw,
            gradients[0].yaw.nodes[node_index][component],
        ),
        (
            AeroOutput::InducedDrag,
            gradients[0].induced_drag.nodes[node_index][component],
        ),
    ] {
        let finite = (mesh_coefficient(&plus, &condition, output)
            - mesh_coefficient(&minus, &condition, output))
            / (2.0 * step);
        assert_abs_diff_eq!(analytical, finite, epsilon = 2.0e-6);
    }
}

#[test]
fn adjoint_batched_transpose_identity() {
    let aircraft = rectangular_wing(2, 4);
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&FlightCondition::default()).unwrap();
    let rhs = VecX::from_element(solver.mesh().vortex_count(), 0.25);
    let lambda = pegasus_vlm::gradients::AdjointSolver::new(&solver)
        .unwrap()
        .solve(&rhs)
        .unwrap();
    assert_abs_diff_eq!(
        solver.linearized_aic().transpose() * lambda,
        rhs,
        epsilon = 1.0e-10
    );
}

#[test]
fn invalid_solution_is_rejected_by_adjoint() {
    let solver = VortexLatticeSolver::new();
    assert!(
        TotalDerivative::compute(&solver, &FlightCondition::default(), AeroOutput::Lift).is_err()
    );
    let state = SolutionState::default();
    assert!(!state.is_valid);
}

/// Perturbs one section design variable on a clone.
fn perturbed(
    aircraft: &pegasus_vlm::geometry::Aircraft,
    surface_index: usize,
    section_index: usize,
    variable: GeometryVariable,
    delta: f64,
) -> pegasus_vlm::geometry::Aircraft {
    let mut clone = aircraft.clone();
    let section = &mut clone.surfaces_mut()[surface_index].sections_mut()[section_index];
    match variable {
        GeometryVariable::Chord => section.set_chord(section.chord() + delta),
        GeometryVariable::Twist => section.set_twist_degrees(section.twist_degrees() + delta),
        GeometryVariable::LeadingEdgeX
        | GeometryVariable::LeadingEdgeY
        | GeometryVariable::LeadingEdgeZ => {
            let mut point = *section.leading_edge();
            let component = match variable {
                GeometryVariable::LeadingEdgeX => 0,
                GeometryVariable::LeadingEdgeY => 1,
                _ => 2,
            };
            point[component] += delta;
            section.set_leading_edge(point);
        }
    }
    clone
}

const GEOMETRY_VARIABLES: [GeometryVariable; 5] = [
    GeometryVariable::Chord,
    GeometryVariable::Twist,
    GeometryVariable::LeadingEdgeX,
    GeometryVariable::LeadingEdgeY,
    GeometryVariable::LeadingEdgeZ,
];

const ALL_OUTPUTS: [AeroOutput; 7] = [
    AeroOutput::Lift,
    AeroOutput::Drag,
    AeroOutput::Side,
    AeroOutput::Roll,
    AeroOutput::Pitch,
    AeroOutput::Yaw,
    AeroOutput::InducedDrag,
];

/// Condition with every rate non-zero, so roll and yaw are genuinely excited.
fn rich_condition() -> FlightCondition {
    FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        roll_rate: 0.03,
        pitch_rate: -0.02,
        yaw_rate: 0.01,
        ..FlightCondition::default()
    }
}

fn check_section_gradients(
    aircraft: &pegasus_vlm::geometry::Aircraft,
    outputs: &[AeroOutput],
    label: &str,
) {
    let condition = rich_condition();
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(aircraft);
    let gradients = calculator.compute_geometry_gradients(&condition).unwrap();
    let section_count = aircraft.surfaces()[0].sections().len();
    assert_eq!(gradients.sections.len(), section_count);

    for section_index in 0..section_count {
        for variable in GEOMETRY_VARIABLES {
            let step = if variable == GeometryVariable::Twist {
                1.0e-4
            } else {
                1.0e-6
            };
            let plus = perturbed(aircraft, 0, section_index, variable, step);
            let minus = perturbed(aircraft, 0, section_index, variable, -step);
            for &output in outputs {
                let analytical = gradients.sections[section_index].get(variable).get(output);
                let finite = (coefficient(&plus, &condition, output)
                    - coefficient(&minus, &condition, output))
                    / (2.0 * step);
                let tolerance = 2.0e-5_f64.max(2.0e-5 * analytical.abs());
                assert!(
                    (analytical - finite).abs() <= tolerance,
                    "{label}: section {section_index} {variable:?} {output:?}: \
                     analytical {analytical:.12} vs finite {finite:.12}"
                );
            }
        }
    }
}

/// The analytical section-geometry chain, verified against full re-solves.
///
/// This is the test that validates replacing the finite-difference producer.
/// The fixture is tapered *and* twisted with a non-zero surface angle offset
/// and a non-unit scale, so a dropped `scale.x`, a dropped `angle_offset`, or
/// the old small-angle twist linearization all fail here.
#[test]
fn section_gradients_match_resolves_on_a_tapered_twisted_wing() {
    check_section_gradients(&common::tapered_twisted_wing(), &ALL_OUTPUTS, "tapered");
}

/// Mirrored geometry, restricted to the outputs that expose sign errors.
///
/// Lift and induced drag are symmetric and would pass even with the mirror
/// reflection wrong; roll and yaw would not.
#[test]
fn section_gradients_match_resolves_on_a_mirrored_wing() {
    check_section_gradients(
        &common::tapered_twisted_wing_mirrored(),
        &[
            AeroOutput::Lift,
            AeroOutput::Roll,
            AeroOutput::Yaw,
            AeroOutput::InducedDrag,
        ],
        "mirrored",
    );
}

/// The only fixture that makes `claf` design-dependent, so the only one that
/// exercises the CLAF-dependent control-point locations.
#[test]
fn section_gradients_handle_a_design_dependent_cl_alpha_factor() {
    for spacing in [
        pegasus_vlm::core::SpacingType::Uniform,
        pegasus_vlm::core::SpacingType::Cosine,
        pegasus_vlm::core::SpacingType::Sine,
        pegasus_vlm::core::SpacingType::Blend,
    ] {
        let mut aircraft = common::cl_alpha_factor_wing();
        aircraft.surfaces_mut()[0].set_chordwise_spacing(spacing);
        for section in aircraft.surfaces_mut()[0].sections_mut() {
            section.set_airfoil(pegasus_vlm::geometry::Airfoil::naca4(2412, 101).unwrap());
        }
        check_section_gradients(
            &aircraft,
            &[AeroOutput::Lift, AeroOutput::Pitch, AeroOutput::InducedDrag],
            "CLAF control-point and camber motion",
        );
    }
}

/// The exact chord-weighted twist derivative, against central differences of
/// the primal `atan2`.
#[test]
fn interpolated_twist_jvp_matches_centered_differences() {
    use pegasus_vlm::gradients::interpolated_twist_jvp;

    let (fc, cl, cr) = (0.37, 1.6, 0.5);
    let (tl, tr) = (6.0f64.to_radians(), (-5.0f64).to_radians());
    let primal = |cl: f64, cr: f64, tl: f64, tr: f64| {
        ((1.0 - fc) * cl * tl.sin() + fc * cr * tr.sin())
            .atan2((1.0 - fc) * cl * tl.cos() + fc * cr * tr.cos())
    };
    let step = 1.0e-7;

    // One seed at a time, in the argument order of the kernel.
    let cases: [(f64, f64, f64, f64); 4] = [
        (1.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0),
        (0.0, 0.0, 0.0, 1.0),
    ];
    for (d_cl, d_cr, d_tl, d_tr) in cases {
        let analytical = interpolated_twist_jvp(fc, cl, cr, tl, tr, d_cl, d_cr, d_tl, d_tr);
        let finite = (primal(
            cl + step * d_cl,
            cr + step * d_cr,
            tl + step * d_tl,
            tr + step * d_tr,
        ) - primal(
            cl - step * d_cl,
            cr - step * d_cr,
            tl - step * d_tl,
            tr - step * d_tr,
        )) / (2.0 * step);
        assert_abs_diff_eq!(analytical, finite, epsilon = 1.0e-8);
    }
}

/// Guards against reintroducing the v0.2.0 small-angle twist linearization.
///
/// The C++ `PropagateTwistGradient` distributed the twist sensitivity as
/// `(1-fc)` and `fc`. That is exact only when the chords *and* twists match at
/// both ends of the interval; on a tapered, twisted interval it is off by tens
/// of percent, which no tolerance would absorb.
#[test]
fn interpolated_twist_jvp_is_not_the_small_angle_linearization() {
    use pegasus_vlm::gradients::interpolated_twist_jvp;

    let fc = 0.5;
    let (cl, cr) = (1.6, 0.5);
    let (tl, tr) = (6.0f64.to_radians(), (-5.0f64).to_radians());
    let exact = interpolated_twist_jvp(fc, cl, cr, tl, tr, 0.0, 0.0, 1.0, 0.0);
    let linearized = 1.0 - fc;
    assert!(
        (exact - linearized).abs() > 0.2,
        "exact {exact} is too close to the linearized {linearized} to discriminate"
    );
    assert_abs_diff_eq!(exact, 0.7638, epsilon = 1.0e-3);
}

/// With equal twists the interpolated angle *is* that twist, so the two
/// endpoint sensitivities must sum to exactly one whatever the chords are.
#[test]
fn interpolated_twist_jvp_sums_to_one_for_uniform_twist() {
    use pegasus_vlm::gradients::interpolated_twist_jvp;

    for fc in [0.0, 0.25, 0.5, 0.75, 1.0] {
        for (cl, cr) in [(1.0, 1.0), (1.6, 0.5), (0.3, 2.2)] {
            for twist in [0.0, 4.0f64.to_radians(), (-7.0f64).to_radians()] {
                let left = interpolated_twist_jvp(fc, cl, cr, twist, twist, 0.0, 0.0, 1.0, 0.0);
                let right = interpolated_twist_jvp(fc, cl, cr, twist, twist, 0.0, 0.0, 0.0, 1.0);
                assert_abs_diff_eq!(left + right, 1.0, epsilon = 1.0e-14);
                if (cl - cr).abs() > 1.0e-12 && fc > 0.0 && fc < 1.0 {
                    // Chord weighting, not linear weighting.
                    assert!((left - (1.0 - fc)).abs() > 1.0e-6);
                }
            }
        }
    }
}

/// The mesher must not move panel vertices when twist changes.
///
/// The whole twist path rests on this: `d(bound_vortex)/d(twist)` and
/// `d(control_point)/d(twist)` are taken as exactly zero. If anyone ever makes
/// the mesher rotate panels about the quarter chord, this fires instead of the
/// gradients silently going wrong.
#[test]
fn twist_perturbation_leaves_panel_vertices_exactly_unchanged() {
    use pegasus_vlm::solver::VortexLattice;

    let aircraft = common::tapered_twisted_wing();
    let mut twisted = aircraft.clone();
    let section = &mut twisted.surfaces_mut()[0].sections_mut()[1];
    section.set_twist_degrees(section.twist_degrees() + 2.0);

    let mut base = VortexLattice::default();
    base.generate_from_aircraft(&aircraft).unwrap();
    let mut moved = VortexLattice::default();
    moved.generate_from_aircraft(&twisted).unwrap();

    assert_eq!(base.vortex_count(), moved.vortex_count());
    let mut normals_changed = 0;
    for (before, after) in base.vortices().iter().zip(moved.vortices()) {
        assert_eq!(before.bound_vortex_left, after.bound_vortex_left);
        assert_eq!(before.bound_vortex_right, after.bound_vortex_right);
        assert_eq!(before.bound_vortex_center, after.bound_vortex_center);
        assert_eq!(before.control_point, after.control_point);
        assert_eq!(before.local_chord, after.local_chord);
        if before.normal_at_control_point != after.normal_at_control_point {
            normals_changed += 1;
        }
    }
    for (before, after) in base.strips().iter().zip(moved.strips()) {
        assert_eq!(before.chord_center, after.chord_center);
    }
    assert!(normals_changed > 0, "twist must still rotate the normals");
}

/// Sideslip and all three rate derivatives, verified against full re-solves.
///
/// These were previously checked only against the frozen oracle, whose only
/// beta entry is `CY_beta = 0.0` on a planar rectangular wing — zero by
/// symmetry no matter what the kernel does. The fixture here is asymmetric
/// (tapered, twisted, dihedralled) so every term is genuinely non-zero.
#[test]
fn beta_and_rate_gradients_match_resolves() {
    let aircraft = common::tapered_twisted_wing();
    let base = rich_condition();
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_flight_condition_gradients(&base)
        .unwrap();

    let step = 1.0e-6;
    let perturb = |field: usize, delta: f64| {
        let mut condition = base.clone();
        match field {
            0 => condition.beta += delta,
            1 => condition.roll_rate += delta,
            2 => condition.pitch_rate += delta,
            _ => condition.yaw_rate += delta,
        }
        condition
    };

    for (field, name, group) in [
        (0usize, "beta", gradients.by_beta),
        (1, "roll rate", gradients.by_roll_rate),
        (2, "pitch rate", gradients.by_pitch_rate),
        (3, "yaw rate", gradients.by_yaw_rate),
    ] {
        for &output in &ALL_OUTPUTS {
            let analytical = group.get(output);
            let finite = (coefficient(&aircraft, &perturb(field, step), output)
                - coefficient(&aircraft, &perturb(field, -step), output))
                / (2.0 * step);
            assert!(
                (analytical - finite).abs() <= 3.0e-5_f64.max(3.0e-5 * analytical.abs()),
                "{name} {output:?}: analytical {analytical:.12} vs finite {finite:.12}"
            );
        }
        // A structurally-zero group would make the comparison vacuous.
        assert!(
            ALL_OUTPUTS.iter().any(|&o| group.get(o).abs() > 1.0e-6),
            "{name} gradients are all zero; the check would be vacuous"
        );
    }
}

/// Section drag polars must appear in the gradients, as they do in the primal.
///
/// `strip_loads` folds a polar-interpolated viscous term into the strip's body
/// force, which then flows into the total force and, through `arm x force`,
/// into the moments. Lift and side are untouched because the added vector is
/// parallel to the stability drag axis and cancels exactly; drag, roll, pitch
/// and yaw are not. That asymmetry is asserted here, so a viscous term applied
/// to the wrong outputs fails just as loudly as a missing one.
#[test]
fn viscous_polars_appear_in_flight_condition_gradients() {
    let aircraft = common::viscous_tapered_wing();
    let base = rich_condition();
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let gradients = calculator
        .compute_flight_condition_gradients(&base)
        .unwrap();

    // The fixture must actually engage the viscous branch.
    let inviscid = common::tapered_twisted_wing();
    assert!(
        (coefficient(&aircraft, &base, AeroOutput::Drag)
            - coefficient(&inviscid, &base, AeroOutput::Drag))
        .abs()
            > 1.0e-3,
        "fixture carries no viscous drag; the test would be vacuous"
    );

    let step = 1.0e-6;
    for (field, name, group) in [
        (0usize, "alpha", gradients.by_alpha),
        (1, "beta", gradients.by_beta),
    ] {
        for &output in &ALL_OUTPUTS {
            let mut plus = base.clone();
            let mut minus = base.clone();
            if field == 0 {
                plus.alpha += step;
                minus.alpha -= step;
            } else {
                plus.beta += step;
                minus.beta -= step;
            }
            let analytical = group.get(output);
            let finite = (coefficient(&aircraft, &plus, output)
                - coefficient(&aircraft, &minus, output))
                / (2.0 * step);
            assert!(
                (analytical - finite).abs() <= 3.0e-5_f64.max(3.0e-5 * analytical.abs()),
                "{name} {output:?}: analytical {analytical:.12} vs finite {finite:.12}"
            );
        }
    }
}

/// Geometry gradients refuse a viscous mesh rather than answering without the
/// polar term.
///
/// The flight-condition path carries the viscous chain, but the geometry path
/// would need strip-area and quarter-chord sensitivities that
/// `LatticeDirections` does not represent. An inviscid-only answer for a
/// viscous mesh is a silently wrong derivative, so it is refused.
#[test]
fn geometry_gradients_reject_a_viscous_mesh() {
    let aircraft = common::viscous_tapered_wing();
    let mut calculator = GradientCalculator::new();
    calculator.set_aircraft(&aircraft);
    let error = calculator
        .compute_geometry_gradients(&rich_condition())
        .unwrap_err();
    assert!(
        matches!(error, pegasus_vlm::PegasusError::InvalidGeometry(_)),
        "expected a typed refusal, got {error:?}"
    );

    // The same wing without polars is still fully supported.
    let mut inviscid = GradientCalculator::new();
    inviscid.set_aircraft(&common::tapered_twisted_wing());
    assert!(
        inviscid
            .compute_geometry_gradients(&rich_condition())
            .is_ok()
    );
}
