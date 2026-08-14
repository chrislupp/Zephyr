//! Ported from the frozen v0.2.0 `mass_reader_test.cpp`.
//!
//! The AVL `.mass` format: unit multipliers, constants, per-component rows with
//! optional products of inertia, and the `*` scale and `+` offset prefixes that
//! apply to every row after them. Aggregation converts to SI on the way out, so
//! the unit multipliers are load-bearing rather than decorative.

use approx::assert_abs_diff_eq;
use pegasus_vlm::io::MassReader;

const TOL: f64 = 1.0e-9;

/// Ports `EmptyContent`.
///
/// Empty input yields no components and unit multipliers of one, so a file that
/// declares nothing behaves as if it were already SI.
#[test]
fn empty_content_gives_unit_multipliers() {
    let data = MassReader::read_str("").unwrap();
    assert_eq!(data.components.len(), 0);
    assert_abs_diff_eq!(data.length_unit, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(data.mass_unit, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(data.time_unit, 1.0, epsilon = TOL);
}

/// Ports `OnlyComments`.
#[test]
fn a_file_of_only_comments_has_no_components() {
    let data = MassReader::read_str(
        "\n\
# This is a comment
! Another comment
  # Indented comment
",
    )
    .unwrap();
    assert_eq!(data.components.len(), 0);
}

/// Ports `UnitDefinitions`.
#[test]
fn unit_definitions_are_parsed() {
    let data = MassReader::read_str(
        "\n\
Lunit = 0.0254 m
Munit = 0.001 kg
Tunit = 1.0 s
",
    )
    .unwrap();
    assert_abs_diff_eq!(data.length_unit, 0.0254, epsilon = TOL);
    assert_abs_diff_eq!(data.mass_unit, 0.001, epsilon = TOL);
    assert_abs_diff_eq!(data.time_unit, 1.0, epsilon = TOL);
}

/// Ports `ConstantDefinitions`.
#[test]
fn gravity_and_density_constants_are_parsed() {
    let data = MassReader::read_str("\ng = 9.81\nrho = 1.225\n").unwrap();
    assert_abs_diff_eq!(data.gravity, 9.81, epsilon = TOL);
    assert_abs_diff_eq!(data.density, 1.225, epsilon = TOL);
}

/// Ports `CaseInsensitiveKeywords`.
#[test]
fn keywords_are_case_insensitive() {
    let data = MassReader::read_str(
        "\n\
LUNIT = 2.0 m
munit = 3.0 kg
TuNiT = 4.0 s
G = 10.0
RHO = 1.5
",
    )
    .unwrap();
    assert_abs_diff_eq!(data.length_unit, 2.0, epsilon = TOL);
    assert_abs_diff_eq!(data.mass_unit, 3.0, epsilon = TOL);
    assert_abs_diff_eq!(data.time_unit, 4.0, epsilon = TOL);
    assert_abs_diff_eq!(data.gravity, 10.0, epsilon = TOL);
    assert_abs_diff_eq!(data.density, 1.5, epsilon = TOL);
}

/// Ports `SingleComponent`.
///
/// The trailing comment is the component's name, so the comment marker both
/// terminates the numeric fields and begins the label.
#[test]
fn a_single_component_row_parses_with_its_trailing_name() {
    let data = MassReader::read_str(
        "\n\
# mass   x     y     z     Ixx   Iyy   Izz
  100.0  0.5   0.0   0.0   10.0  50.0  55.0  ! fuselage
",
    )
    .unwrap();
    assert_eq!(data.components.len(), 1);
    let component = &data.components[0];
    assert_abs_diff_eq!(component.mass, 100.0, epsilon = TOL);
    assert_abs_diff_eq!(component.position.x, 0.5, epsilon = TOL);
    assert_abs_diff_eq!(component.position.y, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(component.position.z, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(component.i_xx, 10.0, epsilon = TOL);
    assert_abs_diff_eq!(component.i_yy, 50.0, epsilon = TOL);
    assert_abs_diff_eq!(component.i_zz, 55.0, epsilon = TOL);
    assert_eq!(component.name, "fuselage");
}

/// Ports `ComponentWithProductsOfInertia`.
///
/// The optional trailing three fields are `Ixy`, `Ixz`, `Iyz` in that order —
/// note `Ixz` precedes `Iyz`, which is easy to transpose.
#[test]
fn optional_products_of_inertia_are_read_in_order() {
    let data = MassReader::read_str(
        "\n  100.0  0.5  0.0  0.0  10.0  50.0  55.0  1.0  2.0  3.0  ! with POI\n",
    )
    .unwrap();
    assert_eq!(data.components.len(), 1);
    let component = &data.components[0];
    assert_abs_diff_eq!(component.i_xy, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(component.i_xz, 2.0, epsilon = TOL);
    assert_abs_diff_eq!(component.i_yz, 3.0, epsilon = TOL);
}

/// Ports `MultipleComponents`.
#[test]
fn multiple_components_keep_their_order_and_names() {
    let data = MassReader::read_str(
        "\n\
Lunit = 1.0 m
Munit = 1.0 kg

  100.0  0.5   0.0   0.0   10.0  50.0  55.0  ! fuselage
   50.0  0.3   2.0   0.1    5.0   1.0   5.0  ! right wing
   50.0  0.3  -2.0   0.1    5.0   1.0   5.0  ! left wing
",
    )
    .unwrap();
    assert_eq!(data.components.len(), 3);
    assert_eq!(data.components[0].name, "fuselage");
    assert_eq!(data.components[1].name, "right wing");
    assert_eq!(data.components[2].name, "left wing");
}

/// Ports `ScaleFactors`, `OffsetFactors` and `ScaleAndOffset`.
///
/// Scale applies before offset: `2 * 1 + 5 = 7`, not `2 * (1 + 5) = 12`. That
/// ordering is the whole content of the combined case.
#[test]
fn scale_applies_before_offset() {
    let scaled =
        MassReader::read_str("\n*  2.0  3.0  4.0\n  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! scaled\n")
            .unwrap();
    assert_eq!(scaled.components.len(), 1);
    assert_abs_diff_eq!(scaled.components[0].position.x, 2.0, epsilon = TOL);
    assert_abs_diff_eq!(scaled.components[0].position.y, 3.0, epsilon = TOL);
    assert_abs_diff_eq!(scaled.components[0].position.z, 4.0, epsilon = TOL);

    let offset = MassReader::read_str(
        "\n+  10.0  20.0  30.0\n  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! offset\n",
    )
    .unwrap();
    assert_abs_diff_eq!(offset.components[0].position.x, 11.0, epsilon = TOL);
    assert_abs_diff_eq!(offset.components[0].position.y, 21.0, epsilon = TOL);
    assert_abs_diff_eq!(offset.components[0].position.z, 31.0, epsilon = TOL);

    let both = MassReader::read_str(
        "\n\
*  2.0  2.0  2.0
+  5.0  5.0  5.0
  1.0  1.0  1.0  1.0  0.0  0.0  0.0  ! scaled and offset
",
    )
    .unwrap();
    assert_abs_diff_eq!(both.components[0].position.x, 7.0, epsilon = TOL);
    assert_abs_diff_eq!(both.components[0].position.y, 7.0, epsilon = TOL);
    assert_abs_diff_eq!(both.components[0].position.z, 7.0, epsilon = TOL);
}

/// Ports `TotalMass` and `TotalMassWithUnits`.
///
/// The total is converted to SI through `Munit`, so 1000 rows of grams weigh
/// one kilogram.
#[test]
fn total_mass_sums_and_converts_to_si() {
    let data = MassReader::read_str(
        "\n\
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp1
   50.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp2
   25.0  0.0  0.0  0.0  0.0  0.0  0.0  ! comp3
",
    )
    .unwrap();
    assert_abs_diff_eq!(data.total_mass(), 175.0, epsilon = TOL);

    let data = MassReader::read_str(
        "\nMunit = 0.001 kg\n  1000.0  0.0  0.0  0.0  0.0  0.0  0.0  ! grams\n",
    )
    .unwrap();
    assert_abs_diff_eq!(data.total_mass(), 1.0, epsilon = TOL);
}

/// Ports `CenterOfGravity` and `CenterOfGravityWeighted`.
///
/// The second case is the one that matters: 100 kg at the origin and 50 kg at
/// `x = 3` puts the centroid at `x = 1`, not at the midpoint 1.5.
#[test]
fn the_center_of_gravity_is_mass_weighted() {
    let data = MassReader::read_str(
        "\n\
Lunit = 1.0 m
Munit = 1.0 kg
  50.0  0.0  0.0  0.0  0.0  0.0  0.0  ! at origin
  50.0  2.0  0.0  0.0  0.0  0.0  0.0  ! at x=2
",
    )
    .unwrap();
    let cg = data.center_of_gravity();
    assert_abs_diff_eq!(cg.x, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(cg.y, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(cg.z, 0.0, epsilon = TOL);

    let data = MassReader::read_str(
        "\n\
Lunit = 1.0 m
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  0.0  0.0  0.0  ! 100 kg at origin
   50.0  3.0  0.0  0.0  0.0  0.0  0.0  ! 50 kg at x=3
",
    )
    .unwrap();
    assert_abs_diff_eq!(data.center_of_gravity().x, 1.0, epsilon = TOL);
}

/// Ports `ComputeTotalProperties`.
#[test]
fn total_properties_of_a_single_component_at_the_origin() {
    let data = MassReader::read_str(
        "\n\
Lunit = 1.0 m
Munit = 1.0 kg
  100.0  0.0  0.0  0.0  10.0  20.0  25.0  ! single component at origin
",
    )
    .unwrap();
    let properties = data.total_properties();
    assert_abs_diff_eq!(properties.mass, 100.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.center_of_gravity.x, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.center_of_gravity.y, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.center_of_gravity.z, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_xx, 10.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_yy, 20.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_zz, 25.0, epsilon = TOL);
}

/// Ports `ParallelAxisTheorem`.
///
/// Two point masses of 50 kg at `y = +/-1` have no inertia of their own, so the
/// totals are entirely transport terms: `Ixx` and `Izz` pick up `m d^2` twice
/// over, and `Iyy` stays zero because the offset is along `y` itself. This is
/// the sharpest check in the file — it fails if the transport axis is confused.
#[test]
fn the_parallel_axis_transport_lands_on_the_right_axes() {
    let data = MassReader::read_str(
        "\n\
Lunit = 1.0 m
Munit = 1.0 kg
  50.0  0.0   1.0  0.0  0.0  0.0  0.0  ! at y=+1
  50.0  0.0  -1.0  0.0  0.0  0.0  0.0  ! at y=-1
",
    )
    .unwrap();
    let properties = data.total_properties();
    assert_abs_diff_eq!(properties.mass, 100.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.center_of_gravity.x, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.center_of_gravity.y, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_xx, 100.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_yy, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(properties.i_zz, 100.0, epsilon = TOL);
}

/// Ports `EmptyComponentsGiveZeroMass`.
///
/// With no mass there is no centroid to divide for, so the accessor must return
/// the origin rather than divide by zero.
#[test]
fn an_empty_component_list_has_zero_mass_and_no_centroid() {
    let data = MassReader::read_str("").unwrap();
    assert_abs_diff_eq!(data.total_mass(), 0.0, epsilon = TOL);
    assert_abs_diff_eq!(data.center_of_gravity().norm(), 0.0, epsilon = TOL);
}

/// Ports `CompleteAVLStyleFile`.
///
/// A realistic file with both comment markers, unit conversions and five
/// components. The total is in grams under `Munit = 0.001`, so 780 g is 0.78 kg.
#[test]
fn a_complete_avl_style_mass_file() {
    let data = MassReader::read_str(
        "\n\
#-------------------------------------------------
#  Dimensional unit and calculation setup
#-------------------------------------------------
Lunit = 0.0254 m
Munit = 0.001 kg
Tunit = 1.0 s

g = 9.81
rho = 1.225

#-------------------------------------------------
#  Mass data
#-------------------------------------------------
#  mass    x       y       z     Ixx     Iyy     Izz
   500.0   10.0    0.0     0.0   100.0   200.0   250.0  ! fuselage
   100.0    5.0   20.0     1.0    10.0    50.0    55.0  ! right wing
   100.0    5.0  -20.0     1.0    10.0    50.0    55.0  ! left wing
    50.0   40.0    0.0     2.0     5.0    10.0    12.0  ! horizontal tail
    30.0   40.0    0.0     5.0     2.0     8.0     8.0  ! vertical tail
",
    )
    .unwrap();
    assert_abs_diff_eq!(data.length_unit, 0.0254, epsilon = TOL);
    assert_abs_diff_eq!(data.mass_unit, 0.001, epsilon = TOL);
    assert_abs_diff_eq!(data.gravity, 9.81, epsilon = TOL);
    assert_abs_diff_eq!(data.density, 1.225, epsilon = TOL);
    assert_eq!(data.components.len(), 5);
    assert_abs_diff_eq!(data.total_mass(), 0.78, epsilon = 1.0e-6);
}
