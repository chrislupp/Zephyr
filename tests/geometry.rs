use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{DEGREES_TO_RADIANS, SpacingType, SymmetryType, Vec3};
use pegasus_vlm::geometry::{
    Aircraft, Airfoil, Body, BodySection, ControlSurfaceDefinition, DesignVariableDefinition,
    MeshInput, Section, Surface, SurfaceMesh,
};
use pegasus_vlm::utils::{
    CubicSpline, SplineBoundaryCondition, blended_spacing, body_to_stability_rotation,
    cosine_spacing, rotate_about_x, rotate_about_y, rotate_about_z, sine_spacing, spacing,
    uniform_spacing, wind_to_body_rotation,
};

#[test]
fn flat_plate_has_zero_camber_and_thickness() {
    let airfoil = Airfoil::flat_plate();
    assert!(airfoil.is_valid());
    assert_eq!(airfoil.name(), "Flat Plate");
    for x in [0.0, 0.2, 0.5, 1.0] {
        assert_eq!(airfoil.camber(x), 0.0);
        assert_eq!(airfoil.thickness_at(x), 0.0);
    }
}

#[test]
fn naca_2412_has_expected_geometry() {
    let airfoil = Airfoil::naca4_str("2412", 101).unwrap();
    assert_eq!(airfoil.name(), "NACA 2412");
    assert_abs_diff_eq!(airfoil.max_camber(), 0.02, epsilon = 1.0e-14);
    assert_abs_diff_eq!(airfoil.max_camber_position(), 0.4, epsilon = 1.0e-14);
    assert_abs_diff_eq!(airfoil.max_thickness(), 0.12, epsilon = 2.0e-3);
    assert!(airfoil.thickness_at(0.3) > 0.05);
}

#[test]
fn naca_validation_rejects_bad_codes_and_point_counts() {
    assert!(Airfoil::naca4_str("", 20).is_err());
    assert!(Airfoil::naca4_str("24A2", 20).is_err());
    assert!(Airfoil::naca4(12_000, 20).is_err());
    assert!(Airfoil::naca4(12, 1).is_err());
}

#[test]
fn explicit_camber_line_is_spline_backed() {
    let airfoil = Airfoil::from_coordinates(&[0.0, 0.5, 1.0], &[0.0, 0.1, 0.0], true).unwrap();
    assert_abs_diff_eq!(airfoil.camber(0.5), 0.1, epsilon = 1.0e-12);
    assert_eq!(airfoil.thickness_distribution(), &[0.0, 0.0, 0.0]);
    assert_eq!(airfoil.camber(-1.0), 0.0);
    assert_abs_diff_eq!(airfoil.camber(2.0), 0.0, epsilon = 1.0e-14);
}

#[test]
fn closed_airfoil_coordinates_extract_camber_and_thickness() {
    let x = [1.0, 0.5, 0.0, 0.5, 1.0];
    let y = [0.0, 0.1, 0.0, -0.1, 0.0];
    let airfoil = Airfoil::from_coordinates(&x, &y, false).unwrap();
    assert_abs_diff_eq!(airfoil.camber(0.5), 0.0, epsilon = 1.0e-12);
    assert_abs_diff_eq!(airfoil.thickness_at(0.5), 0.2, epsilon = 1.0e-6);
}

#[test]
fn section_defaults_and_owned_airfoil_are_consistent() {
    let mut section = Section::new(Vec3::new(1.0, 2.0, 3.0), 2.0, 5.0);
    assert_abs_diff_eq!(section.twist_radians(), 5.0 * DEGREES_TO_RADIANS);
    assert!(!section.has_airfoil());
    assert_eq!(section.airfoil().camber(0.5), 0.0);
    section.set_airfoil(Airfoil::naca4(2412, 30).unwrap());
    assert!(section.has_airfoil());
    assert!(section.airfoil().camber(0.4) > 0.0);
}

#[test]
fn section_control_design_and_viscous_metadata_preserve_order() {
    let mut section = Section::default();
    section.add_control_surface(ControlSurfaceDefinition {
        name: "flap".into(),
        ..ControlSurfaceDefinition::default()
    });
    section.add_control_surface(ControlSurfaceDefinition {
        name: "flap".into(),
        gain: 2.0,
        ..ControlSurfaceDefinition::default()
    });
    section.add_design_variable(DesignVariableDefinition {
        name: "twist".into(),
        weight: 0.5,
    });
    section.set_viscous_polar([-1.0, 0.03, 0.0, 0.01, 1.0, 0.03]);
    assert_eq!(section.control_surfaces().len(), 2);
    assert_eq!(section.control_surfaces()[1].gain, 2.0);
    assert_eq!(section.design_variables()[0].name, "twist");
    assert_eq!(section.viscous_polar()[3], 0.01);
}

fn tapered_surface() -> Surface {
    let mut surface = Surface::new("Wing");
    surface.add_section(Section::new(Vec3::new(0.0, 0.0, 0.0), 2.0, 0.0));
    surface.add_section(Section::new(Vec3::new(0.0, 5.0, 0.0), 1.0, 0.0));
    surface
}

#[test]
fn tapered_surface_area_span_and_mean_chord_match_legacy_formulas() {
    let surface = tapered_surface();
    assert_abs_diff_eq!(surface.planform_area(), 7.5);
    assert_abs_diff_eq!(surface.span(), 5.0);
    assert_abs_diff_eq!(surface.mean_aerodynamic_chord(), 5.0 / 3.0);
}

#[test]
fn y_duplicate_doubles_area_and_span() {
    let mut surface = tapered_surface();
    surface.set_y_duplicate(0.0);
    assert_abs_diff_eq!(surface.planform_area(), 15.0);
    assert_abs_diff_eq!(surface.span(), 10.0);
    surface.clear_y_duplicate();
    assert_abs_diff_eq!(surface.planform_area(), 7.5);
}

#[test]
fn numeric_spacing_conventions_match_avl() {
    let mut surface = Surface::default();
    surface.set_numeric_chordwise_spacing(0.0);
    assert_eq!(surface.chordwise_spacing(), SpacingType::Uniform);
    surface.set_numeric_chordwise_spacing(1.0);
    assert_eq!(surface.chordwise_spacing(), SpacingType::Cosine);
    surface.set_numeric_spanwise_spacing(-1.0);
    assert_eq!(surface.spanwise_spacing(), SpacingType::Sine);
}

#[test]
fn aircraft_aggregates_mesh_counts_and_names() {
    let mut surface = tapered_surface();
    surface.set_chordwise_panel_count(3);
    surface.set_spanwise_panel_count(5);
    surface.set_y_duplicate(0.0);
    for section in surface.sections_mut() {
        section.add_control_surface(ControlSurfaceDefinition {
            name: "aileron".into(),
            ..ControlSurfaceDefinition::default()
        });
        section.add_design_variable(DesignVariableDefinition {
            name: "twist".into(),
            weight: 1.0,
        });
    }
    let mut aircraft = Aircraft::new("Test");
    aircraft.add_surface(surface);
    assert_eq!(aircraft.total_vortex_count(), 30);
    assert_eq!(aircraft.total_strip_count(), 10);
    assert_eq!(aircraft.control_surface_names(), vec!["aileron"]);
    assert_eq!(aircraft.design_variable_names(), vec!["twist"]);
    assert!(aircraft.surface_named("Wing").is_some());
}

#[test]
fn aircraft_reference_symmetry_and_body_properties_round_trip() {
    let mut aircraft = Aircraft::default();
    aircraft.set_reference_area(12.0);
    aircraft.set_reference_chord(1.5);
    aircraft.set_reference_span(9.0);
    aircraft.set_mach_number(0.3);
    aircraft.set_y_symmetry(SymmetryType::Symmetric);
    aircraft.set_z_symmetry(SymmetryType::Antisymmetric);
    aircraft.set_z_symmetry_plane(0.25);
    let mut body = Body::new("Fuselage");
    body.add_section(BodySection {
        center: Vec3::new(1.0, 0.0, 0.0),
        width: 0.4,
        height: 0.5,
    });
    aircraft.add_body(body);
    assert_eq!(aircraft.reference_area(), 12.0);
    assert_eq!(aircraft.mach_number(), 0.3);
    assert_eq!(aircraft.y_symmetry(), SymmetryType::Symmetric);
    assert_eq!(aircraft.bodies()[0].sections().len(), 1);
}

#[test]
fn body_defaults_and_transforms_round_trip() {
    let mut body = Body::new("Pod");
    assert_eq!(body.node_count(), 20);
    assert_eq!(*body.scale(), Vec3::repeat(1.0));
    body.set_node_count(30);
    body.set_node_spacing(-1.0);
    body.set_translation(Vec3::new(1.0, 2.0, 3.0));
    body.set_scale(Vec3::new(2.0, 3.0, 4.0));
    body.set_y_duplicated(true);
    assert_eq!(body.node_count(), 30);
    assert_eq!(body.node_spacing(), -1.0);
    assert_eq!(*body.translation(), Vec3::new(1.0, 2.0, 3.0));
    assert!(body.is_y_duplicated());
}

#[test]
fn surface_mesh_checked_access_and_validation() {
    let mut surface = SurfaceMesh::new("mesh", 2, 3);
    surface.set_node(2, 3, Vec3::new(1.0, 2.0, 3.0)).unwrap();
    assert_eq!(*surface.node(2, 3).unwrap(), Vec3::new(1.0, 2.0, 3.0));
    assert!(surface.set_node(3, 0, Vec3::zeros()).is_err());
    assert!(surface.validate().is_ok());
    surface.nodes.pop();
    assert!(surface.validate().is_err());
}

#[test]
fn mesh_input_rejects_empty_and_nonpositive_references() {
    assert!(MeshInput::default().validate().is_err());
    let mut surface = SurfaceMesh::new("mesh", 1, 1);
    surface.nodes = vec![Vec3::zeros(); 4];
    let input = MeshInput {
        consistent_forces: false,
        surfaces: vec![surface],
        ..MeshInput::default()
    };
    assert!(input.validate().is_ok());
}

#[test]
fn spacing_functions_have_expected_endpoints_and_clustering() {
    for values in [
        uniform_spacing(9),
        cosine_spacing(9),
        sine_spacing(9),
        blended_spacing(9, 0.5),
        spacing(9, SpacingType::Blend),
    ] {
        assert_abs_diff_eq!(values[0], 0.0, epsilon = 1.0e-14);
        assert_abs_diff_eq!(*values.last().unwrap(), 1.0, epsilon = 1.0e-14);
        assert!(values.windows(2).all(|pair| pair[1] >= pair[0]));
    }
    assert!(cosine_spacing(9)[1] < uniform_spacing(9)[1]);
}

#[test]
fn rotations_are_orthogonal_and_follow_right_hand_rule() {
    let vector = Vec3::new(1.0, 0.0, 0.0);
    assert_abs_diff_eq!(
        rotate_about_z(&vector, 90.0 * DEGREES_TO_RADIANS),
        Vec3::y(),
        epsilon = 1.0e-12
    );
    assert_abs_diff_eq!(
        rotate_about_y(&vector, 90.0 * DEGREES_TO_RADIANS),
        -Vec3::z(),
        epsilon = 1.0e-12
    );
    assert_abs_diff_eq!(
        rotate_about_x(&Vec3::y(), 90.0 * DEGREES_TO_RADIANS),
        Vec3::z(),
        epsilon = 1.0e-12
    );
    let wind = wind_to_body_rotation(0.2, -0.1);
    let stability = body_to_stability_rotation(0.2);
    assert_abs_diff_eq!(wind.transpose() * wind, pegasus_vlm::core::Mat3::identity());
    assert_abs_diff_eq!(
        stability.transpose() * stability,
        pegasus_vlm::core::Mat3::identity()
    );
}

#[test]
fn cubic_spline_interpolates_values_and_derivatives() {
    let spline = CubicSpline::build(&[0.0, 0.5, 1.0], &[0.0, 0.25, 1.0]).unwrap();
    assert_abs_diff_eq!(spline.evaluate(0.5).unwrap(), 0.25, epsilon = 1.0e-14);
    let all = spline.evaluate_all(0.3).unwrap();
    assert!(all.0.is_finite() && all.1.is_finite() && all.2.is_finite());
    assert_eq!(spline.len(), 3);
}

#[test]
fn cubic_spline_honors_specified_endpoint_slopes() {
    let spline = CubicSpline::build_with_boundaries(
        &[0.0, 1.0],
        &[0.0, 1.0],
        SplineBoundaryCondition::SpecifiedSlope,
        SplineBoundaryCondition::SpecifiedSlope,
        1.0,
        1.0,
    )
    .unwrap();
    assert_abs_diff_eq!(spline.evaluate_derivative(0.0).unwrap(), 1.0);
    assert_abs_diff_eq!(spline.evaluate_derivative(1.0).unwrap(), 1.0);
    assert!(CubicSpline::build(&[0.0, 0.0], &[0.0, 1.0]).is_err());
}

/// A NaN evaluation parameter must propagate as NaN, not panic.
///
/// `find_interval` returns an index the caller decrements, and NaN satisfies
/// neither of its ordered comparisons, so `partition_point` returned 0 and the
/// subtraction underflowed. The legacy binary search converged to the last
/// interval and returned NaN, which is the behaviour restored here.
#[test]
fn spline_evaluation_at_nan_returns_nan_without_panicking() {
    let spline = CubicSpline::build(&[0.0, 0.5, 1.0], &[0.0, 1.0, 0.0]).unwrap();
    assert!(spline.evaluate(f64::NAN).unwrap().is_nan());
    assert!(spline.evaluate_derivative(f64::NAN).unwrap().is_nan());
    assert!(
        spline
            .evaluate_second_derivative(f64::NAN)
            .unwrap()
            .is_nan()
    );

    // Reachable from the airfoil API, which takes a plain chord fraction.
    let airfoil = Airfoil::naca4_str("2412", 41).unwrap();
    assert!(airfoil.camber(f64::NAN).is_nan());
    assert!(airfoil.camber_slope(f64::NAN).is_nan());

    // Finite parameters are unaffected.
    assert_abs_diff_eq!(spline.evaluate(0.5).unwrap(), 1.0, epsilon = 1.0e-12);
}

/// Geometry that cannot produce a usable lattice must be rejected, not meshed
/// into NaN.
///
/// Each of these previously generated "successfully". Coincident stations gave
/// 2 non-finite control points out of 18, which then poisoned every entry of
/// the AIC.
#[test]
fn degenerate_section_geometry_is_rejected() {
    fn wing(stations: &[f64], chords: &[f64]) -> Aircraft {
        let mut surface = Surface::new("wing");
        for (y, chord) in stations.iter().zip(chords) {
            surface.add_section(Section::new(Vec3::new(0.0, *y, 0.0), *chord, 0.0));
        }
        surface.set_chordwise_panel_count(2);
        surface.set_spanwise_panel_count(9);
        let mut aircraft = Aircraft::new("degenerate");
        aircraft.add_surface(surface);
        aircraft
    }

    let mut solver = pegasus_vlm::solver::VortexLatticeSolver::new();

    // Repeating a spanwise station is a standard AVL spacing-break idiom.
    assert!(
        solver
            .set_aircraft(&wing(&[0.0, 2.0, 2.0, 5.0], &[1.0; 4]))
            .is_err()
    );
    // Zero and negative chord.
    assert!(
        solver
            .set_aircraft(&wing(&[0.0, 5.0], &[0.0, 0.0]))
            .is_err()
    );
    assert!(
        solver
            .set_aircraft(&wing(&[0.0, 5.0], &[1.0, -1.0]))
            .is_err()
    );
    // Sections that double back overlap physically.
    assert!(
        solver
            .set_aircraft(&wing(&[0.0, 5.0, 2.0], &[1.0; 3]))
            .is_err()
    );
    // Non-finite coordinates.
    assert!(
        solver
            .set_aircraft(&wing(&[0.0, f64::NAN], &[1.0; 2]))
            .is_err()
    );

    // A well-formed wing still meshes, and every panel is finite.
    solver
        .set_aircraft(&wing(&[0.0, 2.0, 5.0], &[1.0; 3]))
        .unwrap();
    assert!(
        solver
            .mesh()
            .vortices()
            .iter()
            .all(|vortex| vortex.control_point.iter().all(|v| v.is_finite()))
    );
}

/// The counting helpers must agree with the lattice the mesher builds.
///
/// They multiplied the surface-level panel count by the duplicate factor,
/// ignoring per-section overrides and the per-interval rounding, so they
/// reported 20 strips where the mesher produced 19 and 20 where it produced 7.
#[test]
fn aircraft_counts_match_the_generated_lattice() {
    fn check(aircraft: &Aircraft) {
        let mut lattice = pegasus_vlm::solver::VortexLattice::default();
        lattice.generate_from_aircraft(aircraft).unwrap();
        assert_eq!(aircraft.total_strip_count(), lattice.strip_count());
        assert_eq!(aircraft.total_vortex_count(), lattice.vortex_count());
    }

    // Uneven intervals: 20 panels split over spans of 1 and 2 truncates to
    // 6 + 13 = 19, not 20.
    let mut uneven = Surface::new("uneven");
    for y in [0.0, 1.0, 3.0] {
        uneven.add_section(Section::new(Vec3::new(0.0, y, 0.0), 1.0, 0.0));
    }
    uneven.set_spanwise_panel_count(20);
    uneven.set_chordwise_panel_count(4);
    let mut aircraft = Aircraft::new("uneven");
    aircraft.add_surface(uneven);
    check(&aircraft);

    // A per-section override wins over the surface default outright.
    let mut overridden = Surface::new("overridden");
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0);
    root.set_spanwise_panel_count(7);
    overridden.add_section(root);
    overridden.add_section(Section::new(Vec3::new(0.0, 4.0, 0.0), 1.0, 0.0));
    overridden.set_spanwise_panel_count(20);
    overridden.set_chordwise_panel_count(3);
    let mut aircraft = Aircraft::new("overridden");
    aircraft.add_surface(overridden);
    check(&aircraft);

    // And with mirroring on top.
    aircraft.surfaces_mut()[0].set_y_duplicate(0.0);
    check(&aircraft);
}
