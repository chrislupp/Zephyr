//! Exact homogeneity and flight-direction identities; no numerical derivatives.
use pegasus_vlm::{
    aero::ForceCalculator,
    core::Vec3,
    geometry::{MeshInput, SurfaceMesh},
    gradients::{StripLoadAdjoint, StripLoadCotangent},
    solver::{FlightCondition, VortexLatticeSolver},
};

#[test]
fn velocity_cotangents_obey_quadratic_load_homogeneity_and_angle_chain_rule() {
    let mut surface = SurfaceMesh::new("wing", 3, 4);
    for j in 0..=4 {
        for i in 0..=3 {
            let y = j as f64;
            surface
                .set_node(
                    i,
                    j,
                    Vec3::new(0.1 * y + (1. - 0.05 * y) * i as f64 / 3., y, 0.08 * y),
                )
                .unwrap();
        }
    }
    surface.y_duplicated = true;
    let mesh = MeshInput {
        surfaces: vec![surface],
        ..Default::default()
    };
    let condition = FlightCondition {
        alpha: 0.08,
        beta: 0.02,
        ..Default::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&mesh).unwrap();
    solver.solve(&condition).unwrap();
    let n = solver.mesh().strip_count();
    let points = vec![Vec3::new(0.3, 0., 0.1); n];
    let mut seed = StripLoadCotangent::zeros(n);
    for i in 0..n {
        seed.force[i] = Vec3::new(0.2, 0.1, 1. + i as f64 / 10.);
        seed.moment[i] = Vec3::new(0.1, 0.3, -0.2);
    }
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points).unwrap();
    let velocity = adjoint.velocity_cotangents(&[seed.clone()]).unwrap();
    let angles = adjoint
        .flight_condition_cotangents(&[seed.clone()])
        .unwrap();
    let results = ForceCalculator::new().compute_with_references(
        &mesh.reference,
        solver.mesh(),
        solver.solution(),
        &condition,
    );
    let mut value = 0.;
    for s in &results.surface_loads {
        for strip in &s.strip_loads {
            let f = strip.area * strip.body_force;
            let m = strip.area
                * strip.chord
                * Vec3::new(
                    -strip.body_moment.x,
                    strip.body_moment.y,
                    -strip.body_moment.z,
                )
                + (strip.quarter_chord_position - points[strip.strip_index]).cross(&f);
            value += seed.force[strip.strip_index].dot(&f) + seed.moment[strip.strip_index].dot(&m);
        }
    }
    let direction = VortexLatticeSolver::freestream_direction(&condition);
    let da = Vec3::new(
        -condition.alpha.sin() * condition.beta.cos(),
        0.,
        condition.alpha.cos() * condition.beta.cos(),
    );
    let db = Vec3::new(
        -condition.alpha.cos() * condition.beta.sin(),
        condition.beta.cos(),
        -condition.alpha.sin() * condition.beta.sin(),
    );
    let contraction = |v: &Vec3| {
        velocity[0]
            .force_points
            .iter()
            .chain(&velocity[0].collocation_points)
            .map(|x| x.dot(v))
            .sum::<f64>()
    };
    assert!((contraction(&direction) - 2. * value).abs() < 1e-10 * value.abs());
    assert!((contraction(&da) - angles[0].by_alpha).abs() < 1e-10);
    assert!((contraction(&db) - angles[0].by_beta).abs() < 1e-10);
    assert!(adjoint.velocity_cotangents(&[]).unwrap().is_empty());
    assert!(
        adjoint
            .velocity_cotangents(&[StripLoadCotangent::zeros(n + 1)])
            .is_err()
    );
    // Independent force-point seeds, including non-rigid chordwise work.
    let panel_seed: Vec<_> = (0..solver.mesh().vortex_count())
        .map(|i| Vec3::new(0.1 * (i + 1) as f64, 0.2, 1. / (i + 1) as f64))
        .collect();
    let panel_value: f64 = adjoint
        .panel_forces()
        .unwrap()
        .iter()
        .zip(&panel_seed)
        .map(|(f, b)| f.dot(b))
        .sum();
    let derivatives = adjoint
        .panel_force_derivatives(&mesh, &[panel_seed])
        .unwrap();
    let dilation: f64 = derivatives.mesh[0]
        .surfaces
        .iter()
        .zip(&mesh.surfaces)
        .flat_map(|(bars, surface)| bars.iter().zip(&surface.nodes).map(|(b, x)| b.dot(x)))
        .sum();
    assert!((dilation - 2. * panel_value).abs() < 1e-9 * panel_value.abs());
    let speed: f64 = derivatives.velocity[0]
        .force_points
        .iter()
        .chain(&derivatives.velocity[0].collocation_points)
        .map(|b| b.dot(&direction))
        .sum();
    assert!((speed - 2. * panel_value).abs() < 1e-10 * panel_value.abs());
}
