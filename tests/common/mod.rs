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
