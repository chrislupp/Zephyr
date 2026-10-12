//! Internal trim of a prescribed deformed mesh with exact output sensitivities.
use pegasus_vlm::{
    core::Vec3,
    geometry::{MeshControlSurface, MeshInput, SurfaceMesh},
    solver::FlightCondition,
    trim::{AeroQuantity, MeshTrimProblem, MeshTrimSolver, TrimBalance, TrimVariable},
};

fn main() -> pegasus_vlm::Result<()> {
    let mut mesh = MeshInput {
        consistent_forces: true,
        ..MeshInput::default()
    };
    for (name, xle, span, chord) in [("Wing", 0.0, 3.0, 1.0), ("Tail", 3.0, 1.3, 0.65)] {
        let mut surface = SurfaceMesh::new(name, 2, 3);
        surface.y_duplicated = true;
        for i in 0..=2 {
            for j in 0..=3 {
                let eta = j as f64 / 3.0;
                surface.set_node(
                    i,
                    j,
                    Vec3::new(xle + chord * i as f64 / 2.0, span * eta, 0.08 * eta * eta),
                )?;
            }
        }
        if name == "Tail" {
            surface.controls.push(MeshControlSurface::trailing_edge(
                "elevator",
                0..3,
                vec![1.2; 4],
            ));
        }
        mesh.surfaces.push(surface);
    }
    mesh.reference.area = 7.0;
    mesh.reference.chord = 0.9;
    mesh.reference.span = 6.0;
    mesh.reference.moment_reference_point = Vec3::new(0.3, 0.0, 0.0);
    let mut alpha = TrimVariable::new("alpha", 0.03);
    alpha.bounds = [-0.3, 0.3];
    let mut elevator = TrimVariable::new("control_elevator", 0.0);
    elevator.bounds = [-0.4, 0.4];
    let problem = MeshTrimProblem {
        variables: vec![alpha, elevator],
        balances: vec![
            TrimBalance::new(AeroQuantity::Fz, 1000.0),
            TrimBalance::new(AeroQuantity::My, 1000.0),
        ],
        ..MeshTrimProblem::default()
    };
    let condition = FlightCondition {
        mach: 0.25,
        dynamic_pressure: 650.0,
        ..FlightCondition::default()
    };
    // Prescribed external body-z load; this case enforces Fz and My only.
    let external = [0.0, 0.0, -1200.0, 0.0, 0.0, 0.0];
    let solution = MeshTrimSolver.solve(&mesh, &condition, &problem, &[0.0, 0.0], &external)?;
    let partials = solution.compute_partials(&mesh, &[AeroQuantity::CL, AeroQuantity::CD], true)?;
    println!(
        "alpha = {:.6} rad, elevator = {:.6} rad",
        solution.states[0], solution.states[1]
    );
    println!(
        "CL = {:.6}, CD = {:.6}",
        solution.values[0], solution.values[1]
    );
    println!(
        "physical balance residuals = {:?}",
        solution.diagnostics.residuals
    );
    println!(
        "dCD/d(dynamic_pressure) through trim = {:.9} /Pa",
        partials[&("CD".into(), "dynamic_pressure".into())][0]
    );
    println!(
        "dalpha/d(external_Fz) = {:.9} rad/N",
        partials[&("alpha".into(), "external_Fz".into())][0]
    );
    Ok(())
}
