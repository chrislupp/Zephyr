//! JSON-line bridge from an external optimizer to the explicit Philote server.
//! Run `trim_optimization.py` against a fresh `philote_vlm_server`.
use ndarray::{ArrayD, IxDyn};
use philote_mdo::{ArrayMap, PartialMap, client::ExplicitClient};
use serde_json::{Value, json};
use std::{
    collections::HashMap,
    io::{self, BufRead, Write},
};

const Q: f64 = 650.0;
const TAIL_HALFSPAN: f64 = 1.3;

fn scalar(value: f64) -> ArrayD<f64> {
    ArrayD::from_elem(IxDyn(&[1]), value)
}

fn configuration() -> HashMap<String, Value> {
    HashMap::from([
        ("num_surfaces".into(), json!(2)),
        ("surface_names".into(), json!(["Wing", "Tail"])),
        ("surface_nx".into(), json!([3, 3])),
        ("surface_ny".into(), json!([6, 4])),
        ("surface_y_duplicated".into(), json!([true, true])),
        ("consistent_forces".into(), json!(true)),
        (
            "controls".into(),
            json!([{"surface":"Tail", "name":"elevator", "span_start":0, "span_end":4, "hinge":2.0}]),
        ),
        ("outputs".into(), json!(["CD", "Fz", "My"])),
        (
            "fixed_inputs".into(),
            json!([
                "reference_chord",
                "moment_ref_x",
                "moment_ref_y",
                "moment_ref_z",
                "beta",
                "mach",
                "dynamic_pressure",
                "roll_rate",
                "pitch_rate",
                "yaw_rate",
                "external_Fx",
                "external_Fy",
                "external_Mx",
                "external_My",
                "external_Mz",
                "trim_target_Fz",
                "trim_target_My"
            ]),
        ),
        (
            "trim".into(),
            json!({"enabled":true, "variables":["alpha", "control_elevator"],
            "balances":[{"quantity":"Fz", "scale":1000.0}, {"quantity":"My", "scale":1000.0}],
            "initial_guess":{"alpha":0.04, "control_elevator":-0.06},
            "bounds":{"alpha":[-0.3,0.3], "control_elevator":[-0.5,0.5]},
            "absolute_tolerance":1e-12, "relative_tolerance":0.0}),
        ),
    ])
}

// External sizing relations are differentiated here and composed with the
// component partials. Neither alpha nor elevator is an optimizer design variable.
fn inputs(design: [f64; 2]) -> (ArrayMap, [f64; 2]) {
    let [halfspan, tail_chord] = design;
    let weight = 9.81 * (90.0 + 12.0 * (halfspan / 3.0).powi(3) + 3.0 * tail_chord / 0.65);
    let weight_gradient = [9.81 * 12.0 * halfspan.powi(2) / 9.0, 9.81 * 3.0 / 0.65];
    let mut inputs = HashMap::new();
    for (name, ny, xle, span, chord) in [
        ("Wing", 6, 0.0, halfspan, 1.0),
        ("Tail", 4, 3.0, TAIL_HALFSPAN, tail_chord),
    ] {
        let mut nodes = Vec::new();
        for i in 0..=3 {
            for j in 0..=ny {
                let eta = j as f64 / ny as f64;
                let xi = i as f64 / 3.0;
                nodes.extend([
                    xle + chord * xi,
                    span * eta,
                    0.08 * eta * eta + 0.015 * xi * eta + if name == "Tail" { 0.3 } else { 0.0 },
                ]);
            }
        }
        inputs.insert(
            name.into(),
            ArrayD::from_shape_vec(IxDyn(&[nodes.len()]), nodes).expect("node shape"),
        );
    }
    for (name, value) in [
        ("reference_area", 2.0 * halfspan),
        ("reference_chord", 1.0),
        ("reference_span", 2.0 * halfspan),
        ("moment_ref_x", 0.15),
        ("moment_ref_y", 0.0),
        ("moment_ref_z", 0.0),
        ("beta", 0.0),
        ("mach", 0.25),
        ("dynamic_pressure", Q),
        ("roll_rate", 0.0),
        ("pitch_rate", 0.0),
        ("yaw_rate", 0.0),
        ("external_Fx", 0.0),
        ("external_Fy", 0.0),
        ("external_Fz", -weight),
        ("external_Mx", 0.0),
        ("external_My", 0.0),
        ("external_Mz", 0.0),
        ("trim_target_Fz", 0.0),
        ("trim_target_My", 0.0),
    ] {
        inputs.insert(name.into(), scalar(value));
    }
    (inputs, weight_gradient)
}

fn design_gradient(partials: &PartialMap, output: &str, weight_gradient: [f64; 2]) -> [f64; 2] {
    let block = |input: &str| &partials[&(output.to_string(), input.to_string())];
    let mut gradient = [0.0; 2];
    for i in 0..=3 {
        for j in 0..=6 {
            gradient[0] += block("Wing")[3 * (i * 7 + j) + 1] * j as f64 / 6.0;
        }
    }
    for i in 0..=3 {
        for j in 0..=4 {
            gradient[1] += block("Tail")[3 * (i * 5 + j)] * i as f64 / 3.0;
        }
    }
    gradient[0] += 2.0 * block("reference_area")[0] + 2.0 * block("reference_span")[0];
    for d in 0..2 {
        gradient[d] -= block("external_Fz")[0] * weight_gradient[d];
    }
    gradient
}

async fn evaluate(
    client: &mut ExplicitClient,
    design: [f64; 2],
) -> Result<Value, Box<dyn std::error::Error>> {
    if design.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err("positive finite design variables required".into());
    }
    let (input, weight_gradient) = inputs(design);
    let values = client.compute_function(&input).await?;
    let partials = client.compute_gradient(&input).await?;
    let area = 2.0 * design[0];
    let tail_area = 2.0 * TAIL_HALFSPAN * design[1];
    let cd = values["CD"][0];
    let cd_gradient = design_gradient(&partials, "CD", weight_gradient);
    // Analytical area cost representing profile drag, included in the objective
    // only. The exact inner balance equations remain Fz + external_Fz = 0, My = 0.
    let area_drag = Q * (0.003 * area + 0.006 * tail_area);
    let objective = Q * area * cd + area_drag;
    let objective_gradient = [
        Q * (area * cd_gradient[0] + 2.0 * cd) + Q * 0.006,
        Q * area * cd_gradient[1] + Q * 0.006 * 2.0 * TAIL_HALFSPAN,
    ];
    let elevator = values["control_elevator"][0];
    let elevator_gradient = design_gradient(&partials, "control_elevator", weight_gradient);
    let alpha_gradient = design_gradient(&partials, "alpha", weight_gradient);
    Ok(
        json!({"design":design, "objective":objective, "objective_gradient":objective_gradient,
        "alpha":values["alpha"][0], "alpha_gradient":alpha_gradient,
        "elevator":elevator, "elevator_gradient":elevator_gradient,
        "CD":cd, "area_drag":area_drag, "weight":-input["external_Fz"][0],
        "force_residual":values["Fz"][0] + input["external_Fz"][0], "moment_residual":values["My"][0]}),
    )
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let endpoint = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "http://127.0.0.1:50051".into());
    let mut client = ExplicitClient::connect(endpoint).await?;
    client.set_options(configuration()).await?;
    client.setup().await?;
    client.get_variable_definitions().await?;
    client.get_partial_definitions().await?;
    for line in io::stdin().lock().lines() {
        let design: [f64; 2] = serde_json::from_str(&line?)?;
        let response = match evaluate(&mut client, design).await {
            Ok(value) => value,
            Err(error) => json!({"error":error.to_string()}),
        };
        println!("{response}");
        io::stdout().flush()?;
    }
    Ok(())
}
