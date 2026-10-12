//! Explicit Philote trim registry, derivative, cache, and wire-contract checks.
#![cfg(feature = "philote")]
#[path = "common/mesh_trim.rs"]
mod fixture;
use ndarray::{ArrayD, IxDyn};
use pegasus_vlm::philote::VlmDiscipline;
use pegasus_vlm::trim::AeroQuantity;
use philote_mdo::{
    ArrayMap,
    traits::{Discipline, ExplicitDiscipline},
};
use serde_json::{Value, json};
use std::{collections::HashMap, sync::Arc};

fn array(v: Vec<f64>) -> ArrayD<f64> {
    ArrayD::from_shape_vec(IxDyn(&[v.len()]), v).unwrap()
}

fn configuration() -> HashMap<String, Value> {
    HashMap::from([
        ("num_surfaces".into(), json!(2)),
        ("surface_names".into(), json!(["Wing", "Tail"])),
        ("surface_nx".into(), json!([2, 2])),
        ("surface_ny".into(), json!([3, 3])),
        ("surface_y_duplicated".into(), json!([true, true])),
        ("consistent_forces".into(), json!(true)),
        (
            "controls".into(),
            json!([
                {"surface":"Wing", "name":"flap", "span_start":0, "span_end":3, "hinge":1.2},
                {"surface":"Tail", "name":"elevator", "span_start":0, "span_end":3, "hinge":1.2},
                {"surface":"Tail", "name":"tab", "span_start":1, "span_end":3, "hinge":1.6}
            ]),
        ),
        (
            "outputs".into(),
            json!(AeroQuantity::ALL.map(AeroQuantity::name)),
        ),
        (
            "trim".into(),
            json!({"enabled":true, "variables":["alpha", "control_elevator"],
            "balances":[{"quantity":"Fz", "scale":1000.0}, {"quantity":"My", "scale":1000.0}],
            "initial_guess":{"alpha":0.03, "control_elevator":0.0},
            "bounds":{"alpha":[-0.3, 0.3], "control_elevator":[-0.4, 0.4]},
            "variable_scales":{"alpha":0.1, "control_elevator":0.2}}),
        ),
    ])
}

fn inputs() -> ArrayMap {
    let mesh = fixture::mesh(true);
    let c = fixture::condition();
    let p = fixture::problem();
    let mut inputs = HashMap::new();
    for s in &mesh.surfaces {
        inputs.insert(
            s.name.clone(),
            array(s.nodes.iter().flat_map(|v| [v.x, v.y, v.z]).collect()),
        );
    }
    for (name, value) in [
        ("reference_area", mesh.reference.area),
        ("reference_chord", mesh.reference.chord),
        ("reference_span", mesh.reference.span),
        ("moment_ref_x", mesh.reference.moment_reference_point.x),
        ("moment_ref_y", mesh.reference.moment_reference_point.y),
        ("moment_ref_z", mesh.reference.moment_reference_point.z),
        ("beta", c.beta),
        ("mach", c.mach),
        ("dynamic_pressure", c.dynamic_pressure),
        ("roll_rate", c.roll_rate),
        ("pitch_rate", c.pitch_rate),
        ("yaw_rate", c.yaw_rate),
        ("control_flap", c.control_deflections[1].1),
        ("control_tab", c.control_deflections[2].1),
    ] {
        inputs.insert(name.into(), array(vec![value]));
    }
    for q in &AeroQuantity::ALL[6..] {
        inputs.insert(format!("external_{}", q.name()), array(vec![0.0]));
    }
    for (b, value) in p
        .balances
        .iter()
        .zip(fixture::targets(&mesh, &c, &p, &[0.0; 6]))
    {
        inputs.insert(b.target_input.clone(), array(vec![value]));
    }
    inputs
}

fn discipline(config: &HashMap<String, Value>) -> VlmDiscipline {
    let mut d = VlmDiscipline::default();
    d.initialize().unwrap();
    d.set_options(config).unwrap();
    d.setup().unwrap();
    d.setup_partials().unwrap();
    d
}

#[tokio::test]
async fn trimmed_explicit_outputs_and_partials_match_reconverged_inputs() {
    let d = discipline(&configuration());
    let inputs = inputs();
    // A gradient request must converge trim even when no function request preceded it.
    let partials = d.compute_partials(&inputs).await.unwrap();
    let outputs = d.compute(&inputs).await.unwrap();
    fixture::close(outputs["alpha"][0], fixture::condition().alpha);
    assert_eq!(outputs.len(), 14);
    let meta = d.get_variable_definitions().unwrap();
    for name in ["alpha", "control_elevator"] {
        let v = meta.iter().find(|m| m.name == name).unwrap();
        assert_eq!(
            v.r#type,
            philote_mdo::philote_info::VariableType::KOutput as i32
        );
        assert!(!partials.keys().any(|(_, input)| input == name));
    }
    assert!(
        d.trim_diagnostics()
            .unwrap()
            .unwrap()
            .residuals
            .iter()
            .all(|r| r.abs() < 1e-6)
    );
    assert_eq!(partials, d.compute_partials(&inputs).await.unwrap());
    for (name, index, h) in [
        ("Wing", 17, 1e-5),
        ("Tail", 28, 1e-5),
        ("dynamic_pressure", 0, 1e-3),
        ("trim_target_Fz", 0, 1e-3),
        ("trim_target_My", 0, 1e-3),
        ("external_Fz", 0, 1e-3),
        ("external_My", 0, 1e-3),
        ("control_flap", 0, 1e-5),
        ("moment_ref_x", 0, 1e-5),
    ] {
        let mut plus = inputs.clone();
        let mut minus = inputs.clone();
        plus.get_mut(name).unwrap()[index] += h;
        minus.get_mut(name).unwrap()[index] -= h;
        let p = d.compute(&plus).await.unwrap();
        let m = d.compute(&minus).await.unwrap();
        for output in outputs.keys() {
            fixture::close(
                partials[&(output.clone(), name.into())][index],
                (p[output][0] - m[output][0]) / (2.0 * h),
            );
        }
    }
    let mut invalid = inputs.clone();
    invalid.get_mut("trim_target_Fz").unwrap()[0] = 1e9;
    assert!(d.compute(&invalid).await.is_err());
    assert!(d.trim_diagnostics().unwrap().is_none());
    assert_eq!(outputs, d.compute(&inputs).await.unwrap());
}

#[tokio::test]
async fn coefficient_targets_mixer_and_fixed_inputs_keep_complete_state_derivatives() {
    let mut config = configuration();
    config.insert("outputs".into(), json!(["CD"]));
    let trim = config.get_mut("trim").unwrap();
    trim["variables"] = json!(["alpha", "control_pitch"]);
    trim["balances"] = json!([{"quantity":"CL"}, {"quantity":"Cm"}]);
    trim["initial_guess"] = json!({"alpha":0.03, "control_pitch":0.0});
    trim["bounds"] = json!({"alpha":[-0.3, 0.3], "control_pitch":[-0.4, 0.4]});
    trim["variable_scales"] = json!({"alpha":0.1, "control_pitch":0.2});
    trim["mixer"] = json!({"commands":["control_pitch"], "surfaces":["elevator", "tab"], "matrix":[[1.0], [-0.4]], "offsets":[0.01, 0.02], "bounds":{"elevator":[-0.3, 0.3]}});
    let mut i = inputs();
    i.remove("control_tab");
    i.retain(|n, _| !n.starts_with("external_") && !n.starts_with("trim_target_"));
    i.insert("trim_target_CL".into(), array(vec![0.5]));
    i.insert("trim_target_Cm".into(), array(vec![-0.03]));
    config.insert(
        "fixed_inputs".into(),
        json!(["Wing", "Tail", "mach", "control_flap"]),
    );
    let d = discipline(&config);
    let exact = d.compute_partials(&i).await.unwrap();
    assert!(
        !exact
            .keys()
            .any(|(_, n)| ["Wing", "Tail", "mach", "control_flap"].contains(&n.as_str()))
    );
    assert!(exact.contains_key(&("alpha".into(), "trim_target_CL".into())));
    let outputs = d.compute(&i).await.unwrap();
    assert_eq!(outputs.len(), 5);
    fixture::close(
        outputs["control_tab"][0],
        0.02 - 0.4 * outputs["control_pitch"][0],
    );
    let h = 1e-5;
    let mut p = i.clone();
    let mut m = i.clone();
    p.get_mut("trim_target_CL").unwrap()[0] += h;
    m.get_mut("trim_target_CL").unwrap()[0] -= h;
    let plus = d.compute(&p).await.unwrap();
    let minus = d.compute(&m).await.unwrap();
    for name in outputs.keys() {
        fixture::close(
            exact[&(name.clone(), "trim_target_CL".into())][0],
            (plus[name][0] - minus[name][0]) / (2.0 * h),
        );
    }
}

#[tokio::test]
async fn dimensional_outputs_without_trim_have_physical_units_and_derivatives() {
    let mut config = configuration();
    config.insert("trim".into(), json!({"enabled":false}));
    let d = discipline(&config);
    let mut inputs = inputs();
    inputs.insert("alpha".into(), array(vec![fixture::condition().alpha]));
    inputs.insert("control_elevator".into(), array(vec![0.06]));
    inputs.retain(|n, _| !n.starts_with("external_") && !n.starts_with("trim_target_"));
    let exact = d.compute_partials(&inputs).await.unwrap();
    let values = d.compute(&inputs).await.unwrap();
    assert!(d.trim_diagnostics().unwrap().is_none());
    for q in &AeroQuantity::ALL[6..] {
        fixture::close(
            exact[&(q.name().into(), "dynamic_pressure".into())][0],
            values[q.name()][0] / fixture::condition().dynamic_pressure,
        );
    }
    let mut p = inputs.clone();
    let mut m = inputs.clone();
    p.get_mut("alpha").unwrap()[0] += 1e-6;
    m.get_mut("alpha").unwrap()[0] -= 1e-6;
    let plus = d.compute(&p).await.unwrap();
    let minus = d.compute(&m).await.unwrap();
    for q in AeroQuantity::ALL {
        fixture::close(
            exact[&(q.name().into(), "alpha".into())][0],
            (plus[q.name()][0] - minus[q.name()][0]) / 2e-6,
        );
    }
}

#[tokio::test]
async fn trim_cache_invalidates_force_options_and_rejects_malformed_scalars() {
    let mut d = discipline(&configuration());
    let input = inputs();
    let original = d.compute(&input).await.unwrap();
    d.set_options(&HashMap::from([("consistent_forces".into(), json!(false))]))
        .unwrap();
    assert!(d.trim_diagnostics().unwrap().is_none());
    let legacy = d.compute(&input).await.unwrap();
    fixture::close(legacy["Fz"][0], input["trim_target_Fz"][0]);
    assert!((original["alpha"][0] - legacy["alpha"][0]).abs() > 1e-5);
    let mut invalid = input.clone();
    invalid.insert("mach".into(), array(vec![0.2, 0.3]));
    assert!(d.compute(&invalid).await.is_err());
    invalid.insert("mach".into(), array(vec![f64::NAN]));
    assert!(d.compute_partials(&invalid).await.is_err());
    assert_eq!(legacy, d.compute(&input).await.unwrap());
    d.registry_mut().clear();
    d.setup().unwrap();
    d.setup_partials().unwrap();
    assert_eq!(legacy, d.compute(&input).await.unwrap());
}

#[test]
fn malformed_configuration_and_name_collisions_are_rejected() {
    let bad = [
        json!({"enabled":"yes"}),
        json!({"unexpected":true}),
        json!({"variables":"alpha", "balances":[]}),
        json!({"variables":["alpha"], "balances":[{"quantity":"typo"}]}),
        json!({"variables":["alpha"], "balances":[{"quantity":"CL", "unexpected":0}]}),
        json!({"variables":["alpha"], "balances":[{"quantity":"CL"}], "initial_guess":{"typo":0.1}}),
        json!({"variables":["alpha"], "balances":[{"quantity":"CL"}], "bounds":{"alpha":[0.0]}}),
        json!({"variables":["alpha"], "balances":[{"quantity":"CL"}], "max_iterations":2.5}),
        json!({"variables":["alpha"], "balances":[{"quantity":"CL"}], "mixer":{"commands":["control_pitch"]}}),
    ];
    for trim in bad {
        let mut d = VlmDiscipline::default();
        assert!(
            d.set_options(&HashMap::from([("trim".into(), trim)]))
                .is_err()
        );
    }
    for change in [
        json!({"variables":["alpha", "alpha"]}),
        json!({"variables":["unknown", "control_elevator"]}),
        json!({"balances":[{"quantity":"Fz", "target_input":"alpha"}, {"quantity":"My"}]}),
        json!({"balances":[{"quantity":"Fz", "scale":0.0}, {"quantity":"My"}]}),
        json!({"max_iterations":0}),
        json!({"relative_tolerance":-1.0}),
    ] {
        let mut c = configuration();
        for (k, v) in change.as_object().unwrap() {
            c.get_mut("trim").unwrap()[k] = v.clone();
        }
        let mut d = VlmDiscipline::default();
        assert!(d.set_options(&c).and_then(|()| d.setup()).is_err());
    }
    let mut config = configuration();
    config.insert("fixed_inputs".into(), json!(["alpha"]));
    let mut d = VlmDiscipline::default();
    d.set_options(&config).unwrap();
    assert!(d.setup().is_err());
    let mut d = discipline(&configuration());
    assert!(
        d.set_options(&HashMap::from([("trim".into(), json!({"enabled":false}))]))
            .is_err()
    );
    let mut config = configuration();
    config.insert("trim".into(), Value::String(config["trim"].to_string()));
    assert!(discipline(&config).get_variable_definitions().is_ok());
}

#[tokio::test]
async fn explicit_rpc_discovers_solved_states_and_returns_trimmed_gradients() {
    use philote_mdo::philote_info::{
        discipline_service_server::DisciplineServiceServer,
        explicit_service_server::ExplicitServiceServer,
    };
    use philote_mdo::{client::ExplicitClient, server::ExplicitServer};
    let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
    let address = listener.local_addr().unwrap();
    let server = Arc::new(ExplicitServer::new(VlmDiscipline::default()));
    let (stop, stopped) = tokio::sync::oneshot::channel();
    let task = tokio::spawn(
        tonic::transport::Server::builder()
            .add_service(DisciplineServiceServer::from_arc(server.clone()))
            .add_service(ExplicitServiceServer::from_arc(server))
            .serve_with_incoming_shutdown(
                tokio_stream::wrappers::TcpListenerStream::new(listener),
                async {
                    let _ = stopped.await;
                },
            ),
    );
    let mut client = ExplicitClient::connect(format!("http://{address}"))
        .await
        .unwrap();
    client.set_options(configuration()).await.unwrap();
    client.setup().await.unwrap();
    let meta = client.get_variable_definitions().await.unwrap();
    client.get_partial_definitions().await.unwrap();
    assert!(meta.iter().any(|m| m.name == "alpha"
        && m.r#type == philote_mdo::philote_info::VariableType::KOutput as i32));
    let i = inputs();
    let derivatives = client.compute_gradient(&i).await.unwrap();
    let values = client.compute_function(&i).await.unwrap();
    assert!(derivatives.contains_key(&("alpha".into(), "trim_target_Fz".into())));
    fixture::close(values["alpha"][0], fixture::condition().alpha);
    let mut plus = i.clone();
    plus.get_mut("trim_target_Fz").unwrap()[0] += 0.01;
    let shifted = client.compute_function(&plus).await.unwrap();
    fixture::close(
        (shifted["CD"][0] - values["CD"][0]) / 0.01,
        derivatives[&("CD".into(), "trim_target_Fz".into())][0],
    );
    stop.send(()).unwrap();
    task.await.unwrap().unwrap();
}
