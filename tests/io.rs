use approx::assert_abs_diff_eq;
use pegasus_vlm::core::{SpacingType, SymmetryType, Vec3};
use pegasus_vlm::io::{AvlReader, MassReader};

const SIMPLE_AVL: &str = r#"
Test Wing
0.2
1 -1 0.25
10.0 1.0 10.0
0.5 0.0 0.0
0.012

SURFACE
Wing
4 1.0 8 -1.0
COMPONENT
3
YDUPLICATE
0.0
ANGLE
2.0
SCALE
1.0 2.0 1.0
TRANSLATE
0.5 0.0 0.1
CDCL
-1.0 0.03 0.0 0.01 1.0 0.03

SECTION
0.0 0.0 0.0 1.0 1.0 3 0.0
NACA
2412
CLAF
0.9
CONTROL
elevator 1.5 0.7 0 1 0 -1
DESIGN
twist 0.5

SECTION
0.1 2.0 0.2 0.8 -1.0
NACA
0012
CONTROL
elevator 1.0 0.7 0 1 0 -1
DESIGN
twist 1.0

BODY
Fuselage
20 1.0
SCALE
2.0 1.0 1.0
TRANSLATE
0.1 0.2 0.3
YDUPLICATE
BSEC
0.0 0.0 0.0 0.2 0.3
BSEC
2.0 0.0 0.0 0.5 0.6
"#;

#[test]
fn avl_reader_parses_geometry_modifiers_and_metadata() {
    let aircraft = AvlReader::new().read_str(SIMPLE_AVL).unwrap();
    assert_eq!(aircraft.title(), "Test Wing");
    assert_eq!(aircraft.mach_number(), 0.2);
    assert_eq!(aircraft.y_symmetry(), SymmetryType::Symmetric);
    assert_eq!(aircraft.z_symmetry(), SymmetryType::Antisymmetric);
    assert_eq!(aircraft.z_symmetry_plane(), 0.25);
    assert_eq!(aircraft.profile_drag_coefficient(), 0.012);
    assert_eq!(aircraft.surfaces().len(), 1);
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.chordwise_panel_count(), 4);
    assert_eq!(wing.spanwise_panel_count(), 8);
    assert_eq!(wing.chordwise_spacing(), SpacingType::Cosine);
    assert_eq!(wing.spanwise_spacing(), SpacingType::Sine);
    assert_eq!(wing.component_index(), 3);
    assert_eq!(wing.y_duplicate(), Some(0.0));
    assert_eq!(wing.angle_offset_degrees(), 2.0);
    assert_eq!(*wing.scale(), Vec3::new(1.0, 2.0, 1.0));
    assert_eq!(*wing.translation(), Vec3::new(0.5, 0.0, 0.1));
    assert!(wing.has_viscous_drag());
}

#[test]
fn avl_reader_parses_sections_controls_designs_and_airfoils() {
    let aircraft = AvlReader::new().read_str(SIMPLE_AVL).unwrap();
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.sections().len(), 2);
    assert_eq!(wing.sections()[0].airfoil().name(), "NACA 2412");
    assert_eq!(wing.sections()[1].airfoil().name(), "NACA 12");
    assert_eq!(wing.sections()[0].spanwise_panel_count(), 3);
    assert_eq!(wing.sections()[0].spanwise_spacing(), SpacingType::Uniform);
    assert_eq!(wing.sections()[0].control_surfaces()[0].name, "elevator");
    assert_eq!(wing.sections()[0].control_surfaces()[0].gain, 1.5);
    assert_eq!(wing.sections()[0].design_variables()[0].name, "twist");
    assert_eq!(aircraft.control_surface_names(), vec!["elevator"]);
    assert_eq!(aircraft.design_variable_names(), vec!["twist"]);
}

#[test]
fn avl_reader_parses_body_at_end_of_file() {
    let aircraft = AvlReader::new().read_str(SIMPLE_AVL).unwrap();
    assert_eq!(aircraft.bodies().len(), 1);
    let body = &aircraft.bodies()[0];
    assert_eq!(body.name(), "Fuselage");
    assert_eq!(body.node_count(), 20);
    assert!(body.is_y_duplicated());
    assert_eq!(*body.scale(), Vec3::new(2.0, 1.0, 1.0));
    assert_eq!(body.sections().len(), 2);
    assert_eq!(body.sections()[1].width, 0.5);
}

#[test]
fn avl_reader_loads_repository_fixture() {
    let aircraft = AvlReader::new()
        .read_path("tests/fixtures/rect.avl")
        .unwrap();
    assert_eq!(aircraft.title(), "MACH MDAO AVL");
    assert_eq!(aircraft.surfaces().len(), 1);
    assert_eq!(aircraft.surfaces()[0].sections().len(), 2);
    assert_eq!(aircraft.control_surface_names(), vec!["Elevator"]);
}

#[test]
fn avl_parse_error_carries_line_and_path() {
    let error = AvlReader::new()
        .read_str("Bad\n0\n0 0 0\n1 1 1\n0 0 0\nSURFACE\nWing\n1 0\nSECTION\n0 0\n")
        .unwrap_err()
        .to_string();
    assert!(error.contains("line"));
    assert!(error.contains("SECTION"));

    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("bad.avl");
    std::fs::write(
        &path,
        "Bad\n0\n0 0 0\n1 1 1\n0 0 0\nSURFACE\nWing\n1 0\nSECTION\n0 0\n",
    )
    .unwrap();
    let path_error = AvlReader::new().read_path(&path).unwrap_err().to_string();
    assert!(path_error.contains(path.to_str().unwrap()));
}

#[test]
fn missing_avl_file_returns_typed_io_error() {
    let error = AvlReader::new()
        .read_path("tests/fixtures/does-not-exist.avl")
        .unwrap_err()
        .to_string();
    assert!(error.contains("does-not-exist.avl"));
}

const MASS: &str = r#"
Lunit = 0.5 m
MUNIT = 2.0 kg
tunit = 1.0 s
g = 9.80665
rho = 1.225
* 2.0 3.0 4.0
+ 0.1 0.2 0.3
10 1 2 3 4 5 6 0.1 0.2 0.3 ! wing
20 2 1 0 8 9 10 0.4 0.5 0.6 # fuselage
"#;

#[test]
fn mass_reader_parses_units_scales_offsets_and_products() {
    let data = MassReader::read_str(MASS).unwrap();
    assert_eq!(data.length_unit, 0.5);
    assert_eq!(data.mass_unit, 2.0);
    assert_eq!(data.components.len(), 2);
    assert_eq!(data.components[0].name, "wing");
    assert_eq!(data.components[0].position, Vec3::new(2.1, 6.2, 12.3));
    assert_eq!(data.components[1].i_xz, 0.5);
    assert_abs_diff_eq!(data.total_mass(), 60.0);
}

#[test]
fn mass_reader_combines_center_of_gravity_and_parallel_axis_terms() {
    let data = MassReader::read_str(MASS).unwrap();
    let center = data.center_of_gravity();
    assert!(center.iter().all(|value| value.is_finite()));
    let properties = data.total_properties();
    assert!(properties.is_valid());
    assert_eq!(properties.mass, 60.0);
    assert_abs_diff_eq!(properties.center_of_gravity, center);
    assert!(properties.i_xx > 0.0);
    let inertia = properties.inertia_tensor();
    assert_abs_diff_eq!(inertia, inertia.transpose(), epsilon = 1.0e-12);
}

#[test]
fn empty_mass_file_is_valid_empty_data() {
    let data = MassReader::read_str("# comments only\n").unwrap();
    assert_eq!(data.total_mass(), 0.0);
    assert_eq!(data.center_of_gravity(), Vec3::zeros());
    assert!(!data.total_properties().is_valid());
}

/// A trailing comment must never change how a line is classified.
///
/// The reader used to split on `=` before stripping comments, so any `=`
/// inside a comment captured the line as a unit assignment and discarded the
/// mass component entirely — silently, with no error.
#[test]
fn mass_component_survives_an_equals_sign_inside_a_comment() {
    let with_equals = MassReader::read_str("10 1 2 3 4 5 6 ! main gear (x = 1.2 m)\n").unwrap();
    let without = MassReader::read_str("10 1 2 3 4 5 6 ! main gear\n").unwrap();
    assert_eq!(with_equals.components.len(), 1);
    assert_eq!(without.components.len(), 1);
    assert_abs_diff_eq!(with_equals.total_mass(), 10.0, epsilon = 1.0e-12);
    assert_eq!(with_equals.components[0].position, Vec3::new(1.0, 2.0, 3.0));
}

/// An unrecognized `key = value` line falls through to component parsing
/// rather than consuming the line.
#[test]
fn mass_unknown_key_assignment_does_not_swallow_the_line() {
    let data = MassReader::read_str("Nunit = 3.0\n10 1 2 3 4 5 6\n").unwrap();
    assert_eq!(data.components.len(), 1);
    let known = MassReader::read_str("Lunit = 0.0254\n10 1 2 3 4 5 6\n").unwrap();
    assert_abs_diff_eq!(known.length_unit, 0.0254, epsilon = 1.0e-15);
    assert_eq!(known.components.len(), 1);
}

/// Scale and offset lines must not read numbers out of their own comments,
/// and must leave the previous value untouched when short.
#[test]
fn mass_scale_and_offset_ignore_comment_text() {
    let scaled = MassReader::read_str("* 1.0 1.0 ! 7.0\n10 1 2 3 4 5 6\n").unwrap();
    // Only two numeric fields are present, so the scale stays at unity and z
    // is unchanged rather than picking up the 7.0 from the comment.
    assert_eq!(scaled.components[0].position, Vec3::new(1.0, 2.0, 3.0));

    let full = MassReader::read_str("* 1.0 1.0 2.0 ! note\n10 1 2 3 4 5 6\n").unwrap();
    assert_eq!(full.components[0].position, Vec3::new(1.0, 2.0, 6.0));

    let offset = MassReader::read_str("+ 0.0 0.0 1.0 ! shift\n10 1 2 3 4 5 6\n").unwrap();
    assert_eq!(offset.components[0].position, Vec3::new(1.0, 2.0, 4.0));
}

/// Numeric runs stop at the first token that does not parse, matching the
/// legacy stream extraction, instead of skipping it and pulling later numbers
/// into earlier slots.
#[test]
fn mass_component_stops_at_the_first_non_numeric_field() {
    let data = MassReader::read_str("10 1 2 3 kg 4 5 6\n").unwrap();
    assert!(
        data.components.is_empty(),
        "a line with only four leading numbers is not a complete component"
    );
}

/// A malformed panel count must be a typed parse error, never a panic.
///
/// `value as usize` saturates, so `1e30` became `usize::MAX` and overflowed
/// the mesher's panel arithmetic.
#[test]
fn avl_rejects_out_of_range_panel_counts() {
    let header = "rect\n0.0\n0 0 0\n1.0 1.0 1.0\n0.0 0.0 0.0\n";
    for paneling in ["1e30 1.0 4 1.0", "4 1.0 1e30 1.0", "-3 1.0 4 1.0"] {
        let source = format!("{header}SURFACE\nWing\n{paneling}\nSECTION\n0 0 0 1 0\n");
        let error = AvlReader::new().read_str(&source).unwrap_err();
        assert!(
            matches!(error, pegasus_vlm::error::PegasusError::Parse { .. }),
            "expected a parse error for paneling {paneling:?}, got {error:?}"
        );
    }

    let section = format!(
        "{header}SURFACE\nWing\n4 1.0 4 1.0\nSECTION\n0 -1 0 1 0 1e30 1.0\nSECTION\n0 1 0 1 0\n"
    );
    assert!(AvlReader::new().read_str(&section).is_err());
}

/// AVL leaves `SgnDup` off the CONTROL line routinely; requiring it discarded
/// the whole aircraft.
#[test]
fn avl_accepts_a_control_line_without_the_sign_field() {
    let source = "rect\n0.0\n0 0 0\n1.0 1.0 1.0\n0.0 0.0 0.0\n\
                  SURFACE\nWing\n4 1.0 4 1.0\n\
                  SECTION\n0 0 0 1 0\nCONTROL\nelevator 1.0 0.6 0. 1. 0.\n";
    let aircraft = AvlReader::new().read_str(source).unwrap();
    let controls = aircraft.surfaces()[0].sections()[0].control_surfaces();
    assert_eq!(controls.len(), 1);
    assert_eq!(controls[0].name, "elevator");
    assert_abs_diff_eq!(controls[0].sign, 1.0, epsilon = 1.0e-15);

    // An explicit sign still wins.
    let signed = source.replace("0. 1. 0.\n", "0. 1. 0. -1.0\n");
    let aircraft = AvlReader::new().read_str(&signed).unwrap();
    assert_abs_diff_eq!(
        aircraft.surfaces()[0].sections()[0].control_surfaces()[0].sign,
        -1.0,
        epsilon = 1.0e-15
    );

    // Five fields is still malformed.
    let short = source.replace("elevator 1.0 0.6 0. 1. 0.", "elevator 1.0 0.6 0. 1.");
    assert!(AvlReader::new().read_str(&short).is_err());
}

/// A missing airfoil file must not degrade silently to a flat plate: the solve
/// would run on the wrong camber and report success.
#[test]
fn avl_reports_an_unloadable_airfoil_file() {
    let source = "rect\n0.0\n0 0 0\n1.0 1.0 1.0\n0.0 0.0 0.0\n\
                  SURFACE\nWing\n4 1.0 4 1.0\n\
                  SECTION\n0 0 0 1 0\nAFILE\nnot_a_real_airfoil.dat\n";
    let error = AvlReader::new().read_str(source).unwrap_err();
    assert!(
        matches!(error, pegasus_vlm::error::PegasusError::Parse { .. }),
        "expected a parse error, got {error:?}"
    );
}
