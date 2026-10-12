//! Explicit trim client. First run `philote_vlm_server`, then this example.
use ndarray::{ArrayD, IxDyn};
use pegasus_vlm::core::Vec3;
use pegasus_vlm::geometry::SurfaceMesh;
use philote_mdo::client::ExplicitClient;
use serde_json::json;
use std::collections::HashMap;

fn scalar(value: f64) -> ArrayD<f64> {
    ArrayD::from_elem(IxDyn(&[1]), value)
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let endpoint = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "http://127.0.0.1:50051".into());
    let mut client = ExplicitClient::connect(endpoint).await?;
    client.set_options(HashMap::from([
        ("num_surfaces".into(), json!(2)), ("surface_names".into(), json!(["Wing", "Tail"])),
        ("surface_nx".into(), json!([2, 2])), ("surface_ny".into(), json!([3, 3])),
        ("surface_y_duplicated".into(), json!([true, true])), ("consistent_forces".into(), json!(true)),
        ("controls".into(), json!([{"surface":"Tail", "name":"elevator", "span_start":0, "span_end":3, "hinge":1.2}])),
        ("outputs".into(), json!(["CL", "CD", "Fz", "My"])),
        ("trim".into(), json!({"enabled":true, "variables":["alpha", "control_elevator"],
            "balances":[{"quantity":"Fz", "scale":1000.0}, {"quantity":"My", "scale":1000.0}],
            "initial_guess":{"alpha":0.03, "control_elevator":0.0},
            "bounds":{"alpha":[-0.3,0.3], "control_elevator":[-0.4,0.4]}})),
    ])).await?;
    client.setup().await?;
    client.get_variable_definitions().await?;
    client.get_partial_definitions().await?;
    let mut inputs = HashMap::new();
    for (name, xle, span, chord) in [("Wing", 0.0, 3.0, 1.0), ("Tail", 3.0, 1.3, 0.65)] {
        let mut surface = SurfaceMesh::new(name, 2, 3);
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
        let nodes = surface.nodes.iter().flat_map(|n| [n.x, n.y, n.z]).collect();
        inputs.insert(name.into(), ArrayD::from_shape_vec(IxDyn(&[36]), nodes)?);
    }
    for (name, value) in [
        ("reference_area", 7.0),
        ("reference_chord", 0.9),
        ("reference_span", 6.0),
        ("moment_ref_x", 0.3),
        ("moment_ref_y", 0.0),
        ("moment_ref_z", 0.0),
        ("beta", 0.0),
        ("mach", 0.25),
        ("dynamic_pressure", 650.0),
        ("roll_rate", 0.0),
        ("pitch_rate", 0.0),
        ("yaw_rate", 0.0),
        ("external_Fx", 0.0),
        ("external_Fy", 0.0),
        ("external_Fz", -1200.0),
        ("external_Mx", 0.0),
        ("external_My", 0.0),
        ("external_Mz", 0.0),
        ("trim_target_Fz", 0.0),
        ("trim_target_My", 0.0),
    ] {
        inputs.insert(name.into(), scalar(value));
    }
    let values = client.compute_function(&inputs).await?;
    let gradients = client.compute_gradient(&inputs).await?;
    println!(
        "alpha = {:.6}, elevator = {:.6}",
        values["alpha"][0], values["control_elevator"][0]
    );
    println!(
        "dCD/d(external_Fz) = {:.9}",
        gradients[&("CD".into(), "external_Fz".into())][0]
    );
    inputs.insert("external_Fz".into(), scalar(-1250.0));
    let changed = client.compute_function(&inputs).await?;
    println!(
        "alpha after external-load change = {:.6}",
        changed["alpha"][0]
    );
    Ok(())
}
