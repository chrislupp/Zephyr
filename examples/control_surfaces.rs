//! Controls on an externally deformed mesh, with coefficient and load adjoints.
use pegasus_vlm::aero::ForceCalculator;
use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{MeshControlSurface, MeshInput, SurfaceMesh};
use pegasus_vlm::gradients::{GradientCalculator, StripLoadAdjoint, StripLoadCotangent};
use pegasus_vlm::solver::{FlightCondition, VortexLatticeSolver};

fn main() -> pegasus_vlm::Result<()> {
    let mut wing = SurfaceMesh::new("Wing", 4, 6);
    wing.y_duplicated = true;
    for i in 0..=4 {
        for j in 0..=6 {
            let eta = j as f64 / 6.0;
            let x = i as f64 / 4.0;
            // Representative incoming bending/twist deformation, before actuation.
            wing.set_node(
                i,
                j,
                Vec3::new(x, 3.0 * eta, 0.1 * eta * eta + 0.03 * eta * x),
            )?;
        }
    }
    wing.controls.push(MeshControlSurface::trailing_edge(
        "flap",
        1..6,
        vec![3.0; 7],
    ));
    let mut mesh = MeshInput {
        consistent_forces: true,
        ..MeshInput::default()
    };
    mesh.surfaces.push(wing);
    mesh.reference.area = 6.0;
    mesh.reference.span = 6.0;
    let condition = FlightCondition {
        alpha: 4.0_f64.to_radians(),
        control_deflections: vec![("flap".into(), 5.0_f64.to_radians())],
        ..FlightCondition::default()
    };
    let mut solver = VortexLatticeSolver::new();
    solver.set_mesh(&mesh)?;
    let solution = solver.solve(&condition)?;
    let loads = ForceCalculator::new().compute_with_references(
        &mesh.reference,
        solver.mesh(),
        &solution,
        &condition,
    );
    println!(
        "CL = {:.6}, CD = {:.6}",
        loads.coefficients.lift, loads.coefficients.drag
    );

    let mut gradients = GradientCalculator::new();
    gradients.set_mesh(&mesh)?;
    for control in gradients.compute_control_gradients(&condition)? {
        println!(
            "dCL/d{} = {:.6} per radian",
            control.name, control.by_deflection.coefficients.lift
        );
    }

    // An external coupled solver supplies arbitrary strip-force/moment seeds.
    // Here the functional is the summed body-z force before dynamic pressure.
    let count = solver.mesh().strip_count();
    let points = vec![Vec3::zeros(); count];
    let mut seed = StripLoadCotangent::zeros(count);
    seed.force.fill(Vec3::z());
    let adjoint = StripLoadAdjoint::new(&solver, &condition, &points)?;
    println!(
        "Strip-load control derivatives: {:?}",
        adjoint.control_cotangents(&[seed.clone()])?
    );
    let mesh_gradients = adjoint.mesh_node_cotangents(&mesh, &[seed])?;
    println!(
        "Mesh gradient has {} input nodes",
        mesh_gradients[0].surfaces[0].len()
    );
    Ok(())
}
