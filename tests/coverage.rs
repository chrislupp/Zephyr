mod common;

use pegasus_vlm::aero::{
    AerodynamicResults, BodyAxisDerivatives, ForceCalculator, StabilityDerivativeCalculator,
};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, SpacingType, SymmetryType, Vec3, VecX};
use pegasus_vlm::geometry::{
    Airfoil, MeshInput, ReferenceQuantities, Section, Surface, SurfaceMesh,
};
use pegasus_vlm::gradients::{
    AdjointSolver, AeroOutput, ForceCoefficientGradients, GradientCalculator, OutputSelection,
    aic_column_gradient, aic_row_gradient, force_geometry_jvp, force_induced_velocity_gradient,
    selected_force_circulation_gradients, velocity_column_gradient, velocity_row_gradient,
};
use pegasus_vlm::solver::{
    FlightCondition, InfluenceMatrixBuilder, SymmetryConfig, VortexLatticeSolver,
};
use pegasus_vlm::trim::{ConstraintType, TrimPresets};
use pegasus_vlm::utils::{
    CubicSpline, SplineBoundaryCondition, cosine_spacing, prandtl_glauert_factor, sine_spacing,
    trailing_vortex_velocity, uniform_spacing, vortex_segment_velocity,
};
use pegasus_vlm::{PegasusError, Result};

fn solved_wing() -> (
    pegasus_vlm::geometry::Aircraft,
    FlightCondition,
    VortexLatticeSolver,
) {
    let aircraft = common::rectangular_wing(2, 4);
    let condition = FlightCondition {
        alpha: 4.0 * DEGREES_TO_RADIANS,
        beta: 1.5 * DEGREES_TO_RADIANS,
        mach: 0.2,
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    solver.solve(&condition).unwrap();
    (aircraft, condition, solver)
}

#[test]
fn analytical_row_column_and_force_paths_are_finite() {
    let (_, condition, solver) = solved_wing();
    let count = solver.mesh().vortex_count();
    let weights = VecX::from_iterator(count, (0..count).map(|i| i as f64 + 1.0));
    let (control, normal) = aic_row_gradient(solver.mesh(), 0, &weights, 0.98, 0.01, 0.1);
    let (left, right) = aic_column_gradient(solver.mesh(), 0, &weights, 0.98, 0.01, 0.1);
    assert!(
        control
            .iter()
            .chain(normal.iter())
            .all(|value| value.is_finite())
    );
    assert!(
        left.iter()
            .chain(right.iter())
            .all(|value| value.is_finite())
    );

    let cotangent = Vec3::new(0.3, -0.2, 0.8);
    let row = velocity_row_gradient(
        solver.mesh(),
        0,
        &solver.solution().circulation,
        &cotangent,
        0.98,
        0.01,
        0.1,
    );
    let cotangents = vec![cotangent; count];
    let (left, right) =
        velocity_column_gradient(solver.mesh(), 0, 1.0, &cotangents, 0.98, 0.01, 0.1);
    assert!(
        row.iter()
            .chain(left.iter())
            .chain(right.iter())
            .all(|x| x.is_finite())
    );
    assert_eq!(
        velocity_column_gradient(solver.mesh(), 0, 0.0, &cotangents, 0.98, 0.01, 0.1),
        (Vec3::zeros(), Vec3::zeros())
    );

    let selected =
        selected_force_circulation_gradients(&solver, &condition, OutputSelection::all());
    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
    ] {
        assert_eq!(selected.get(output).unwrap().len(), count);
    }
    assert!(selected.get(AeroOutput::InducedDrag).is_none());
    let empty = ForceCoefficientGradients::default();
    assert_eq!(empty.lift.len(), 0);

    for output in [
        AeroOutput::Lift,
        AeroOutput::Drag,
        AeroOutput::Side,
        AeroOutput::Roll,
        AeroOutput::Pitch,
        AeroOutput::Yaw,
        AeroOutput::InducedDrag,
    ] {
        let velocity = force_induced_velocity_gradient(
            solver.mesh(),
            solver.solution(),
            &condition,
            solver.reference_quantities(),
            output,
        );
        assert_eq!(velocity.len(), count);
    }
    let directions = vec![Vec3::new(0.1, -0.2, 0.3); count];
    let geometry = force_geometry_jvp(
        solver.mesh(),
        solver.solution(),
        &condition,
        solver.reference_quantities(),
        &directions,
        &directions,
        &directions,
    )
    .unwrap();
    assert!(geometry.iter().all(|value| value.is_finite()));

    // A slice that does not match the lattice is an error, not a silently
    // truncated partial derivative.
    let short = vec![Vec3::zeros(); count - 1];
    assert!(
        force_geometry_jvp(
            solver.mesh(),
            solver.solution(),
            &condition,
            solver.reference_quantities(),
            &short,
            &directions,
            &directions,
        )
        .is_err()
    );

    let invalid = pegasus_vlm::solver::SolutionState::default();
    assert!(
        force_induced_velocity_gradient(
            solver.mesh(),
            &invalid,
            &condition,
            solver.reference_quantities(),
            AeroOutput::Lift,
        )
        .iter()
        .all(|gradient| *gradient == Vec3::zeros())
    );
}

#[test]
fn influence_builder_exercises_both_symmetry_planes_and_velocity_queries() {
    let (_, _, solver) = solved_wing();
    let symmetry = SymmetryConfig {
        y_symmetry: SymmetryType::Symmetric,
        y_symmetry_plane: 0.0,
        z_symmetry: SymmetryType::Antisymmetric,
        z_symmetry_plane: -0.25,
    };
    let mut builder = InfluenceMatrixBuilder::new(solver.mesh());
    builder.set_mach_number(0.3);
    builder.set_symmetry(symmetry);
    builder.set_core_radius_chord_fraction(0.02);
    builder.set_core_radius_width_fraction(0.03);
    assert_eq!(builder.mach_number(), 0.3);
    assert!(builder.prandtl_glauert_factor() < 1.0);
    assert_eq!(builder.symmetry().z_symmetry, SymmetryType::Antisymmetric);
    let velocity = builder.velocity_influence_matrix();
    let normalwash = builder.normalwash_matrix();
    assert_eq!(
        velocity.shape(),
        (
            3 * solver.mesh().vortex_count(),
            solver.mesh().vortex_count()
        )
    );
    assert_eq!(
        normalwash.shape(),
        (solver.mesh().vortex_count(), solver.mesh().vortex_count())
    );
    let point = Vec3::new(1.5, 0.4, 0.7);
    assert!(builder.single_vortex_velocity(&point, 0, true).norm() > 0.0);
    assert!(
        builder
            .induced_velocity(&point, &solver.solution().circulation)
            .iter()
            .all(|value| value.is_finite())
    );

    builder.set_mach_number(1.2);
    assert_eq!(builder.prandtl_glauert_factor(), 1.0);
}

#[test]
fn body_axis_stability_and_result_helpers_are_populated() {
    let (aircraft, condition, mut solver) = solved_wing();
    let derivative = StabilityDerivativeCalculator
        .compute_body_axis(&aircraft, &mut solver, &condition)
        .unwrap();
    assert!(derivative.is_valid);
    assert!(
        derivative
            .by_alpha_force
            .iter()
            .all(|value| value.is_finite())
    );
    let defaults = BodyAxisDerivatives::default();
    assert!(!defaults.is_valid);
    assert_eq!(defaults.reference_area, 1.0);

    let result = AerodynamicResults {
        reference_area: 2.0,
        reference_span: 4.0,
        ..AerodynamicResults::default()
    };
    assert_eq!(result.aspect_ratio(), 8.0);
    assert_eq!(result.ideal_induced_drag(), 0.0);

    let mut force = ForceCalculator::new();
    force.set_trefftz_core_radius(0.02, 0.2);
    let computed =
        force.compute_for_aircraft(&aircraft, solver.mesh(), solver.solution(), &condition);
    assert!(computed.is_valid);
}

#[test]
fn utility_geometry_and_error_branches_have_typed_behavior() -> Result<()> {
    assert!(CubicSpline::default().evaluate(0.0).is_err());
    assert!(CubicSpline::build(&[0.0], &[0.0]).is_err());
    assert!(CubicSpline::build(&[0.0, 1.0], &[0.0]).is_err());
    assert!(CubicSpline::build(&[0.0, 0.0], &[0.0, 1.0]).is_err());
    let spline = CubicSpline::build_with_boundaries(
        &[0.0, 1.0],
        &[0.0, 1.0],
        SplineBoundaryCondition::ZeroThirdDerivative,
        SplineBoundaryCondition::ZeroThirdDerivative,
        0.0,
        0.0,
    )?;
    assert!(spline.is_valid());
    assert_eq!(spline.len(), 2);
    assert!(!spline.is_empty());
    assert_eq!(spline.min_parameter(), 0.0);
    assert_eq!(spline.max_parameter(), 1.0);
    assert!(spline.evaluate_derivative(0.5)?.is_finite());
    assert!(spline.evaluate_second_derivative(0.5)?.is_finite());

    assert_eq!(uniform_spacing(0), Vec::<f64>::new());
    assert_eq!(uniform_spacing(1), vec![0.0]);
    assert_eq!(cosine_spacing(1), vec![0.0]);
    assert_eq!(sine_spacing(1), vec![0.0]);
    assert_eq!(prandtl_glauert_factor(0.0), 1.0);
    let point = Vec3::new(0.5, 0.5, 0.5);
    let start = Vec3::zeros();
    let end = Vec3::x();
    assert!(vortex_segment_velocity(&point, &start, &end, 0.01).norm() > 0.0);
    assert_eq!(
        vortex_segment_velocity(&start, &start, &end, 0.01),
        Vec3::zeros()
    );
    assert!(trailing_vortex_velocity(&point, &start, 0.01).norm() > 0.0);
    assert_eq!(
        trailing_vortex_velocity(&start, &start, 0.01),
        Vec3::zeros()
    );

    let directory = tempfile::tempdir().unwrap();
    let coordinate_path = directory.path().join("foil.dat");
    std::fs::write(
        &coordinate_path,
        "Coverage foil\n1.0 0.0\n0.5 0.1\n0.0 0.0\n0.5 -0.1\n1.0 0.0\n",
    )
    .unwrap();
    let mut foil = Airfoil::from_coordinate_file(&coordinate_path)?;
    assert_eq!(foil.name(), "Coverage foil");
    foil.set_name("renamed");
    assert_eq!(foil.name(), "renamed");
    assert_eq!(foil.camber_x().len(), foil.camber_y().len());
    assert!(Airfoil::from_coordinate_file(directory.path().join("missing.dat")).is_err());
    let bad_path = directory.path().join("bad.dat");
    std::fs::write(&bad_path, "bad\n0 0\n1 0\n").unwrap();
    assert!(Airfoil::from_coordinate_file(&bad_path).is_err());

    let io_error = PegasusError::io(
        "missing",
        std::io::Error::new(std::io::ErrorKind::NotFound, "not found"),
    );
    assert!(io_error.to_string().contains("missing"));
    assert!(PegasusError::parse(7, "bad").to_string().contains("line 7"));
    assert!(
        PegasusError::parse_path(&coordinate_path, 9, "bad")
            .to_string()
            .contains("line 9")
    );
    Ok(())
}

#[test]
fn geometry_accessors_validation_and_gradient_configuration_are_checked() {
    let mut surface = Surface::new("original");
    surface.set_name("renamed");
    surface.set_component_index(4);
    surface.set_chordwise_spacing(SpacingType::Uniform);
    surface.set_spanwise_spacing(SpacingType::Sine);
    surface.set_angle_offset_degrees(2.0);
    surface.set_translation(Vec3::new(1.0, 2.0, 3.0));
    surface.set_scale(Vec3::new(2.0, 1.0, 0.5));
    surface.set_has_wake(false);
    surface.set_contributes_to_loads(false);
    surface.set_affected_by_alpha_beta(false);
    surface.set_cl_cd_polar([-1.0, 0.1, 0.0, 0.01, 1.0, 0.1]);
    surface.add_section(Section::default());
    assert_eq!(surface.name(), "renamed");
    assert_eq!(surface.component_index(), 4);
    assert_eq!(surface.section(0).unwrap().chord(), 1.0);
    assert!(surface.section_mut(0).is_some());
    assert!(!surface.has_wake());
    assert!(!surface.contributes_to_loads());
    assert!(!surface.affected_by_alpha_beta());
    assert!(surface.has_viscous_drag());
    assert_eq!(surface.span(), 0.0);

    let zero = SurfaceMesh::new("zero", 0, 1);
    assert!(zero.validate().is_err());
    let mut nonfinite = SurfaceMesh::new("bad", 1, 1);
    nonfinite.nodes[0].x = f64::NAN;
    assert!(nonfinite.validate().is_err());
    let invalid_reference = MeshInput {
        surfaces: vec![SurfaceMesh::new("mesh", 1, 1)],
        reference: ReferenceQuantities {
            area: 0.0,
            ..ReferenceQuantities::default()
        },
    };
    assert!(invalid_reference.validate().is_err());

    let mut gradients = GradientCalculator::new();
    assert!(
        gradients
            .set_reference_point(Vec3::new(1.0, 0.0, 0.0))
            .is_err()
    );
    assert!(gradients.set_mesh(&MeshInput::default()).is_err());
    let mesh = common::rectangular_mesh(1, 2);
    gradients.set_mesh(&mesh).unwrap();
    gradients
        .set_reference_point(Vec3::new(0.3, 0.0, 0.0))
        .unwrap();
    gradients.set_output_selection(OutputSelection::only(AeroOutput::Side));
    let result = gradients
        .compute_beta_gradients(&FlightCondition::default())
        .unwrap();
    assert!(result.coefficients.side.is_finite());

    let mut aircraft_gradients = GradientCalculator::new();
    let aircraft = common::rectangular_wing(1, 2);
    aircraft_gradients.set_aircraft(&aircraft);
    aircraft_gradients
        .set_reference_point(Vec3::new(0.2, 0.0, 0.0))
        .unwrap();
    let alpha = aircraft_gradients
        .compute_alpha_gradients(&FlightCondition::default())
        .unwrap();
    assert!(alpha.coefficients.lift.is_finite());

    let level = TrimPresets::level_trimmed_flight(0.3, "elevator");
    assert_eq!(level.controls[0].1.kind, ConstraintType::PitchMoment);
    let turn = TrimPresets::steady_turn(0.4, 0.2);
    assert_eq!(turn.alpha.kind, ConstraintType::Lift);

    let (_, _, solver) = solved_wing();
    let adjoint = AdjointSolver::new(&solver).unwrap();
    assert_eq!(adjoint.unknown_count(), solver.mesh().vortex_count());
    assert!(adjoint.solve(&VecX::zeros(1)).is_err());
    assert!(
        adjoint
            .solve_batched(&pegasus_vlm::core::MatX::zeros(1, 2))
            .is_err()
    );
}
