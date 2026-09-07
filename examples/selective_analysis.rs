//! Run with `cargo run --release --example selective_analysis`.
use pegasus_vlm::analysis::{
    Analysis, AnalysisRequest, GradientRequest, Input, Output, SolverMode,
};
use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::{MeshInput, SurfaceMesh};
use pegasus_vlm::solver::FlightCondition;

fn main() -> pegasus_vlm::Result<()> {
    let mut wing = SurfaceMesh::new("wing", 4, 16);
    wing.y_duplicated = true;
    for i in 0..=wing.nx {
        for j in 0..=wing.ny {
            wing.set_node(i, j, Vec3::new(i as f64 / 4.0, 5.0 * j as f64 / 16.0, 0.0))?;
        }
    }
    let mut input = MeshInput::default();
    input.surfaces.push(wing);
    input.reference.area = 10.0;
    input.reference.span = 10.0;
    let condition = FlightCondition {
        alpha: 5.0_f64.to_radians(),
        ..FlightCondition::default()
    };
    let request = AnalysisRequest {
        mode: SolverMode::HalfWing,
        outputs: vec![Output::Lift, Output::FarFieldInducedDrag],
        gradients: vec![GradientRequest {
            output: Output::Lift,
            input: Input::Alpha,
        }],
        ..AnalysisRequest::default()
    };
    let mut analysis = Analysis::from_mesh(&input);
    let result = analysis.evaluate(&condition, &request)?;
    println!(
        "CL = {:.10}, CDi = {:.10}",
        result.outputs[&Output::Lift],
        result.outputs[&Output::FarFieldInducedDrag]
    );
    println!("dCL/dalpha: {:?}", result.gradients[0].value);
    println!("Work: {:?}", result.work);
    println!(
        "Repeated request: {:?}",
        analysis.evaluate(&condition, &request)?.work
    );
    Ok(())
}
