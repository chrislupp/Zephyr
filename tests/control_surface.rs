//! Ported from the frozen v0.2.0 `control_surface_test.cpp`.
//!
//! Control surfaces travel from the AVL file through `Section`, `Aircraft` and
//! the lattice, ending as a per-panel influence with a hinge axis and a gain.
//! These cases follow that chain and assert what survives each hop.

use approx::assert_abs_diff_eq;
use pegasus_vlm::aero::{ControlDerivatives, StabilityDerivatives};
use pegasus_vlm::core::{DEGREES_TO_RADIANS, Vec3};
use pegasus_vlm::geometry::{Aircraft, ControlSurfaceDefinition, Section, Surface};
use pegasus_vlm::io::AvlReader;
use pegasus_vlm::solver::{FlightCondition, VortexLattice, VortexLatticeSolver};

const TOL: f64 = 1.0e-10;

fn control(
    name: &str,
    hinge_x_fraction: f64,
    hinge_vector: Vec3,
    sign: f64,
) -> ControlSurfaceDefinition {
    ControlSurfaceDefinition {
        name: name.into(),
        hinge_x_fraction,
        hinge_vector,
        gain: 1.0,
        sign,
    }
}

fn aileron() -> ControlSurfaceDefinition {
    control("aileron", 0.75, Vec3::y(), 1.0)
}

/// A wing carrying one named control on both of its sections.
fn wing_with(
    name: &str,
    chordwise: usize,
    spanwise: usize,
    tip_y: f64,
    duplicated: bool,
    sign: f64,
) -> Aircraft {
    let definition = control(name, 0.75, Vec3::y(), sign);
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(chordwise);
    wing.set_spanwise_panel_count(spanwise);
    if duplicated {
        wing.set_y_duplicate(0.0);
    }
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.0, 0.0);
    root.add_control_surface(definition.clone());
    let mut tip = Section::new(Vec3::new(0.0, tip_y, 0.0), 1.0, 0.0);
    tip.add_control_surface(definition);
    wing.add_section(root);
    wing.add_section(tip);

    let mut aircraft = Aircraft::new("Test Aircraft");
    aircraft.add_surface(wing);
    aircraft
}

// ============================================================================
// ControlSurfaceDefinition and Section
// ============================================================================

/// Ports `NameDefaultsToEmpty` and `SetValues`.
#[test]
fn control_surface_definition_fields() {
    assert!(ControlSurfaceDefinition::default().name.is_empty());
    let definition = aileron();
    assert_eq!(definition.name, "aileron");
    assert_abs_diff_eq!(definition.hinge_x_fraction, 0.75, epsilon = TOL);
    assert_abs_diff_eq!(definition.hinge_vector, Vec3::new(0.0, 1.0, 0.0));
    assert_abs_diff_eq!(definition.sign, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(definition.gain, 1.0, epsilon = TOL);
}

/// Ports `NoControlSurfacesByDefault`.
#[test]
fn a_section_has_no_controls_by_default() {
    let section = Section::new(Vec3::zeros(), 1.0, 0.0);
    assert!(section.control_surfaces().is_empty());
}

/// Ports `AddSingleControlSurface`.
#[test]
fn a_section_records_one_control() {
    let mut section = Section::new(Vec3::zeros(), 1.0, 0.0);
    section.add_control_surface(aileron());
    assert_eq!(section.control_surfaces().len(), 1);
    assert_eq!(section.control_surfaces()[0].name, "aileron");
    assert_abs_diff_eq!(
        section.control_surfaces()[0].hinge_x_fraction,
        0.75,
        epsilon = TOL
    );
}

/// Ports `AddMultipleControlSurfaces`.
///
/// Order is preserved, which matters because deflections are applied in list
/// order and an antisymmetric control composes differently than a symmetric
/// one.
#[test]
fn a_section_keeps_multiple_controls_in_order() {
    let mut section = Section::new(Vec3::zeros(), 1.0, 0.0);
    section.add_control_surface(control("flap", 0.70, Vec3::y(), 1.0));
    section.add_control_surface(control("aileron", 0.80, Vec3::y(), -1.0));
    assert_eq!(section.control_surfaces().len(), 2);
    assert_eq!(section.control_surfaces()[0].name, "flap");
    assert_eq!(section.control_surfaces()[1].name, "aileron");
}

// ============================================================================
// Aircraft aggregation
// ============================================================================

/// Ports `GetControlSurfaceNames` and `GetNumControlSurfaces`.
#[test]
fn aircraft_aggregates_one_control_across_sections() {
    let aircraft = wing_with("elevator", 4, 8, 4.0, false, 1.0);
    let names = aircraft.control_surface_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "elevator");
}

/// Ports `MultipleControlSurfaces`.
#[test]
fn aircraft_aggregates_controls_across_surfaces() {
    let mut aircraft = wing_with("elevator", 4, 8, 4.0, false, 1.0);
    let rudder = control("rudder", 0.6, Vec3::z(), 1.0);
    let mut htail = Surface::new("HTail");
    htail.set_chordwise_panel_count(4);
    htail.set_spanwise_panel_count(4);
    let mut root = Section::new(Vec3::new(3.0, 0.0, 0.0), 0.5, 0.0);
    root.add_control_surface(rudder.clone());
    let mut tip = Section::new(Vec3::new(3.0, 1.0, 0.0), 0.5, 0.0);
    tip.add_control_surface(rudder);
    htail.add_section(root);
    htail.add_section(tip);
    aircraft.add_surface(htail);

    let names = aircraft.control_surface_names();
    assert_eq!(names.len(), 2);
    assert!(names.iter().any(|name| name == "elevator"));
    assert!(names.iter().any(|name| name == "rudder"));
}

// ============================================================================
// Lattice propagation
// ============================================================================

fn lattice_for(aircraft: &Aircraft) -> VortexLattice {
    let mut mesh = VortexLattice::new();
    mesh.generate_from_aircraft(aircraft).unwrap();
    mesh
}

/// Ports `ControlSurfaceNamesInMesh` and `ControlInfluencesHaveCorrectSize`.
#[test]
fn control_names_and_influence_rows_reach_the_lattice() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, false, 1.0));
    assert_eq!(mesh.control_surface_names().len(), 1);
    assert_eq!(mesh.control_surface_names()[0], "aileron");
    // One influence list per panel, not per strip or per surface.
    assert_eq!(mesh.control_influences().len(), mesh.vortex_count());
}

/// Ports `TrailingEdgePanelsHaveInfluence`.
#[test]
fn some_panels_carry_a_meaningful_control_gain() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, false, 1.0));
    let influenced = mesh
        .control_influences()
        .iter()
        .flatten()
        .filter(|influence| influence.name == "aileron" && influence.gain.abs() > 0.01)
        .count();
    assert!(influenced > 0);
}

/// Ports `LeadingEdgePanelsNoInfluence`.
///
/// The hinge sits at 75% chord, so the first panel of each strip lies well
/// ahead of it and must not pick up a significant gain. This is what
/// distinguishes a real hinge-line calculation from one that flags the whole
/// strip.
#[test]
fn leading_edge_panels_stay_ahead_of_the_hinge() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, false, 1.0));
    let influences = mesh.control_influences();
    for strip in mesh.strips() {
        for influence in &influences[strip.first_vortex_index] {
            if influence.name == "aileron" {
                assert!(
                    influence.gain.abs() < 0.5,
                    "leading edge panel has gain {}",
                    influence.gain
                );
            }
        }
    }
}

/// Ports `InfluenceHingeAxisCorrect`.
///
/// A hinge vector of `(0, 1, 0)` must survive normalization as a essentially
/// pure spanwise axis.
#[test]
fn the_hinge_axis_stays_spanwise() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, false, 1.0));
    for influence in mesh.control_influences().iter().flatten() {
        if influence.name == "aileron" && influence.gain.abs() > 0.01 {
            assert_abs_diff_eq!(influence.hinge_axis.x, 0.0, epsilon = 0.1);
            assert!(influence.hinge_axis.y.abs() > 0.9);
            assert_abs_diff_eq!(influence.hinge_axis.z, 0.0, epsilon = 0.1);
        }
    }
}

/// Ports `MirroredSurfaceCreated` and `MirroredInfluencesExist`.
#[test]
fn duplication_creates_a_mirrored_surface_with_its_own_influences() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, true, -1.0));
    assert_eq!(mesh.surface_info().len(), 2);
    assert_eq!(mesh.control_influences().len(), mesh.vortex_count());
}

/// Ports `MirroredControlSurfaceInfluences`.
///
/// The mirrored half has to carry the control too. An antisymmetric aileron
/// keeps its magnitude across the mirror even though its sign flips, so the
/// check is on `|gain|`.
#[test]
fn the_mirrored_half_carries_the_same_control() {
    let mesh = lattice_for(&wing_with("aileron", 4, 4, 2.0, true, -1.0));
    let influences = mesh.control_influences();
    let original = &mesh.surface_info()[0];
    let mirror = &mesh.surface_info()[1];

    let mut matched = false;
    for offset in 0..original.num_vortices.min(mirror.num_vortices) {
        let original_row = &influences[original.first_vortex_index + offset];
        let mirror_row = &influences[mirror.first_vortex_index + offset];
        for (left, right) in original_row.iter().zip(mirror_row) {
            if left.name == "aileron" && right.name == "aileron" {
                matched = true;
                if left.gain.abs() > 0.01 {
                    assert!(
                        right.gain.abs() > 0.01,
                        "mirrored panel lost its control influence"
                    );
                }
            }
        }
    }
    assert!(matched, "no aileron influence paired across the mirror");
}

// ============================================================================
// AVL parsing
// ============================================================================

fn avl_with_controls(surface: &str, body: &str) -> String {
    format!(
        "\
Test Aircraft
0.0
0 0 0
10.0 1.0 10.0
0.25 0.0 0.0

SURFACE
{surface}
4 1.0
{body}"
    )
}

/// Ports `ParseBasicControlSurface`.
#[test]
fn avl_parses_a_single_control() {
    let content = avl_with_controls(
        "Wing",
        "\n\
SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
elevator 1.0 0.75 0 1 0 1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
elevator 1.0 0.75 0 1 0 1
",
    );
    let aircraft = AvlReader::default().read_str(&content).unwrap();
    let names = aircraft.control_surface_names();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "elevator");
    let wing = &aircraft.surfaces()[0];
    assert_eq!(wing.sections()[0].control_surfaces().len(), 1);
    assert_eq!(wing.sections()[0].control_surfaces()[0].name, "elevator");
}

/// Ports `ParseMultipleControlSurfaces`.
#[test]
fn avl_parses_two_controls_on_one_section() {
    let content = avl_with_controls(
        "Wing",
        "\n\
SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
flap 1.0 0.70 0 1 0 1
CONTROL
aileron 1.0 0.80 0 1 0 -1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
flap 1.0 0.70 0 1 0 1
CONTROL
aileron 1.0 0.80 0 1 0 -1
",
    );
    let aircraft = AvlReader::default().read_str(&content).unwrap();
    assert_eq!(aircraft.control_surface_names().len(), 2);
    assert_eq!(
        aircraft.surfaces()[0].sections()[0]
            .control_surfaces()
            .len(),
        2
    );
}

/// Ports `ParseControlSurfaceHingeVector`.
///
/// The hinge vector is columns 4-6 of the `CONTROL` line; a rudder's is
/// `(0, 0, 1)`, not the wing default.
#[test]
fn avl_parses_a_vertical_hinge_vector() {
    let content = avl_with_controls(
        "Fin",
        "\n\
SECTION
0.0 0.0 0.0 0.5 0.0
CONTROL
rudder 1.0 0.60 0 0 1 1

SECTION
0.0 0.0 1.0 0.5 0.0
CONTROL
rudder 1.0 0.60 0 0 1 1
",
    );
    let aircraft = AvlReader::default().read_str(&content).unwrap();
    let controls = aircraft.surfaces()[0].sections()[0].control_surfaces();
    assert_eq!(controls.len(), 1);
    assert_eq!(controls[0].name, "rudder");
    assert_abs_diff_eq!(
        controls[0].hinge_vector,
        Vec3::new(0.0, 0.0, 1.0),
        epsilon = TOL
    );
}

/// Ports `ParseControlSurfaceGain`.
///
/// The gain field is signed: a spoiler with `-0.5` deflects opposite its
/// commanded input and at half authority.
#[test]
fn avl_parses_a_negative_control_gain() {
    let content = avl_with_controls(
        "Wing",
        "\n\
SECTION
0.0 0.0 0.0 1.0 0.0
CONTROL
spoiler -0.5 0.50 0 1 0 1

SECTION
0.0 2.0 0.0 1.0 0.0
CONTROL
spoiler -0.5 0.50 0 1 0 1
",
    );
    let aircraft = AvlReader::default().read_str(&content).unwrap();
    let controls = aircraft.surfaces()[0].sections()[0].control_surfaces();
    assert_eq!(controls.len(), 1);
    assert_eq!(controls[0].name, "spoiler");
    assert_abs_diff_eq!(controls[0].gain, -0.5, epsilon = TOL);
}

// ============================================================================
// Solver integration
// ============================================================================

/// Ports `ControlSurfacesPropagateToMesh` and `InfluenceMatrixBuiltWithControls`.
#[test]
fn a_controlled_aircraft_meshes_and_solves() {
    let mut aircraft = wing_with("elevator", 8, 12, 5.0, true, 1.0);
    aircraft.set_reference_area(10.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(10.0);

    let mut solver = VortexLatticeSolver::new();
    solver.set_aircraft(&aircraft).unwrap();
    assert_eq!(solver.mesh().control_surface_names().len(), 1);
    assert_eq!(solver.mesh().control_surface_names()[0], "elevator");

    let solution = solver
        .solve(&FlightCondition {
            alpha: 5.0 * DEGREES_TO_RADIANS,
            ..FlightCondition::default()
        })
        .unwrap();
    assert!(solution.is_valid);
    assert!(!solution.circulation.is_empty());
}

// ============================================================================
// Result types
// ============================================================================

/// Ports `ControlDerivativesTest.DefaultValues` and `SetValues`.
///
/// `gain` defaults to one, not zero: it is a multiplier, and the neutral
/// element for a multiplier is unity. A derived `Default` would silently scale
/// every derivative to nothing.
#[test]
fn control_derivatives_default_to_a_unit_gain() {
    let derivs = ControlDerivatives::default();
    assert!(derivs.name.is_empty());
    assert_abs_diff_eq!(derivs.deflection, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.gain, 1.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.lift, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.side, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.drag, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.roll, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.pitch, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.yaw, 0.0, epsilon = TOL);
    assert_abs_diff_eq!(derivs.hinge_moment, 0.0, epsilon = TOL);

    let derivs = ControlDerivatives {
        name: "elevator".into(),
        coefficients: pegasus_vlm::aero::Coefficients {
            lift: 0.5,
            pitch: -1.2,
            ..pegasus_vlm::aero::Coefficients::default()
        },
        ..ControlDerivatives::default()
    };
    assert_eq!(derivs.name, "elevator");
    assert_abs_diff_eq!(derivs.coefficients.lift, 0.5, epsilon = TOL);
    assert_abs_diff_eq!(derivs.coefficients.pitch, -1.2, epsilon = TOL);
}

/// Ports `FlightConditionControlTest.DefaultNoDeflections` and
/// `SetControlDeflection`.
#[test]
fn flight_condition_control_deflections() {
    assert!(FlightCondition::default().control_deflections.is_empty());
    let condition = FlightCondition {
        control_deflections: vec![
            ("elevator".into(), 5.0 * DEGREES_TO_RADIANS),
            ("aileron".into(), 10.0 * DEGREES_TO_RADIANS),
        ],
        ..FlightCondition::default()
    };
    assert_eq!(condition.control_deflections.len(), 2);
    assert_eq!(condition.control_deflections[0].0, "elevator");
    assert_abs_diff_eq!(
        condition.control_deflections[0].1,
        5.0 * DEGREES_TO_RADIANS,
        epsilon = TOL
    );
}

/// Ports `StabilityDerivativesControlTest.ControlDerivativesVector`.
#[test]
fn stability_derivatives_carry_a_control_list() {
    let mut derivs = StabilityDerivatives::default();
    assert!(derivs.control_derivatives.is_empty());
    derivs.control_derivatives.push(ControlDerivatives {
        name: "elevator".into(),
        ..ControlDerivatives::default()
    });
    assert_eq!(derivs.control_derivatives.len(), 1);
    assert_eq!(derivs.control_derivatives[0].name, "elevator");
}
