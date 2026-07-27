#![allow(dead_code)]

use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{
    Aircraft, ControlSurfaceDefinition, MeshInput, ReferenceQuantities, Section, Surface,
    SurfaceMesh,
};

pub fn rectangular_wing(nx: usize, ny: usize) -> Aircraft {
    let mut wing = Surface::new("Wing");
    wing.set_chordwise_panel_count(nx);
    wing.set_spanwise_panel_count(ny);
    wing.add_section(Section::new(Vec3::new(0.0, -2.0, 0.0), 1.0, 0.0));
    wing.add_section(Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0));
    let mut aircraft = Aircraft::new("Rectangular wing");
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft
}

pub fn controlled_wing() -> Aircraft {
    let control = ControlSurfaceDefinition {
        name: "elevator".into(),
        hinge_x_fraction: 0.7,
        hinge_vector: Vec3::y(),
        sign: 1.0,
        gain: 1.0,
    };
    let mut root = Section::new(Vec3::new(0.0, -2.0, 0.0), 1.0, 0.0);
    root.add_control_surface(control.clone());
    let mut tip = Section::new(Vec3::new(0.0, 2.0, 0.0), 1.0, 0.0);
    tip.add_control_surface(control);
    let mut wing = Surface::new("Controlled wing");
    wing.set_chordwise_panel_count(4);
    wing.set_spanwise_panel_count(6);
    wing.add_section(root);
    wing.add_section(tip);
    let mut aircraft = Aircraft::new("Control fixture");
    aircraft.set_reference_area(4.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(4.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft
}

pub fn rectangular_mesh(nx: usize, ny: usize) -> MeshInput {
    let mut surface = SurfaceMesh::new("Wing", nx, ny);
    for i in 0..=nx {
        for j in 0..=ny {
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(i as f64 / nx as f64, -2.0 + 4.0 * j as f64 / ny as f64, 0.0),
                )
                .unwrap();
        }
    }
    MeshInput {
        surfaces: vec![surface],
        reference: ReferenceQuantities {
            area: 4.0,
            chord: 1.0,
            span: 4.0,
            moment_reference_point: Vec3::new(0.25, 0.0, 0.0),
        },
    }
}

/// Tapered, twisted, swept, dihedralled three-section wing.
///
/// Every interval has `chord_left != chord_right` **and**
/// `twist_left != twist_right`, which is the only regime where the exact
/// chord-weighted twist interpolation differs from the small-angle
/// linearization the v0.2.0 C++ used. A rectangular untwisted wing cannot tell
/// the two apart.
///
/// The explicit per-section spanwise counts are load-bearing for finite-
/// difference verification: without them `per_interval` is
/// `((count * arc_fraction) as usize).max(1)`, an integer truncation that a
/// leading-edge y or z step can cross, silently changing the panel count and
/// producing a garbage oracle.
///
/// A non-zero `angle_offset_degrees` and a non-unit `scale` are set because
/// both are zero/unity in every other fixture, so dropping either from the
/// twist or chord chain would otherwise be invisible.
pub fn tapered_twisted_wing() -> Aircraft {
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.6, 6.0);
    root.set_spanwise_panel_count(3);
    let mut mid = Section::new(Vec3::new(0.30, 1.2, 0.10), 1.1, 1.0);
    mid.set_spanwise_panel_count(3);
    let tip = Section::new(Vec3::new(0.90, 2.4, 0.28), 0.5, -5.0);

    let mut wing = Surface::new("Tapered wing");
    wing.set_chordwise_panel_count(3);
    wing.set_angle_offset_degrees(3.0);
    wing.set_scale(Vec3::new(1.3, 1.0, 1.0));
    wing.add_section(root);
    wing.add_section(mid);
    wing.add_section(tip);

    let mut aircraft = Aircraft::new("Tapered twisted wing");
    aircraft.set_reference_area(3.0);
    aircraft.set_reference_chord(1.0);
    aircraft.set_reference_span(2.4);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft
}

/// [`tapered_twisted_wing`] mirrored about `y = 0`.
///
/// Mirroring sign errors are invisible in lift and induced drag, which are
/// symmetric; they show up in roll and yaw.
pub fn tapered_twisted_wing_mirrored() -> Aircraft {
    let mut aircraft = tapered_twisted_wing();
    aircraft.surfaces_mut()[0].set_y_duplicate(0.0);
    aircraft.set_reference_area(6.0);
    aircraft.set_reference_span(4.8);
    aircraft
}

/// Wing whose sections disagree on `cl_alpha_factor`.
///
/// This is the only configuration that makes `claf` design-dependent and so
/// the only one that exercises `d(x_vortex[0])`. The factors are chosen to keep
/// `0.25 / claf` well inside the `[0.1, 0.4]` clamp, away from its corners.
pub fn cl_alpha_factor_wing() -> Aircraft {
    let mut root = Section::new(Vec3::new(0.0, 0.0, 0.0), 1.4, 2.0);
    root.set_spanwise_panel_count(3);
    root.set_cl_alpha_factor(0.85);
    let mut tip = Section::new(Vec3::new(0.2, 2.0, 0.0), 0.9, -1.0);
    tip.set_cl_alpha_factor(1.0);

    let mut wing = Surface::new("claf wing");
    wing.set_chordwise_panel_count(3);
    wing.add_section(root);
    wing.add_section(tip);

    let mut aircraft = Aircraft::new("cl-alpha factor wing");
    aircraft.set_reference_area(2.3);
    aircraft.set_reference_chord(1.15);
    aircraft.set_reference_span(2.0);
    aircraft.set_moment_reference_point(Vec3::new(0.25, 0.0, 0.0));
    aircraft.add_surface(wing);
    aircraft
}

/// Tapered, twisted wing carrying a section drag polar on every section.
///
/// The polar `cd = 0.01 + 0.02*cl^2` has well-separated breakpoints, so
/// `interpolate_drag_polar` stays on its smooth quadratic branch and away from
/// both the degenerate-abscissa early return and the `max(0.0)` kink.
///
/// No other fixture in the suite sets a viscous polar on a gradient case,
/// which is why the missing viscous term in the gradient module was invisible.
pub fn viscous_tapered_wing() -> Aircraft {
    let mut aircraft = tapered_twisted_wing();
    for section in aircraft.surfaces_mut()[0].sections_mut() {
        section.set_viscous_polar([-1.0, 0.03, 0.0, 0.01, 1.0, 0.03]);
    }
    aircraft
}
