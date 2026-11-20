/*
    Project Zephyr

    Copyright 2024-2026, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-2198.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/

//! Serves the mesh discipline over gRPC.
//!
//! ```text
//! cargo run --features philote --example philote_mesh_server
//! ```
//!
//! Configure the panel counts through discipline options before `setup`:
//!
//! ```json
//! {"nx": 4, "ny": 8, "cosine_chordwise_spacing": true, "y_duplicated": true}
//! ```

use std::sync::Arc;

use philote_mdo::philote_info::{
    discipline_service_server::DisciplineServiceServer,
    explicit_service_server::ExplicitServiceServer,
};
use philote_mdo::server::ExplicitServer;
use tonic::transport::Server;
use zephyr_vlm::philote::MeshDiscipline;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let address = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "127.0.0.1:50052".to_string())
        .parse()?;
    let server = Arc::new(ExplicitServer::new(MeshDiscipline::default()).with_verbose(true));

    println!("Zephyr mesh discipline listening on {address}");

    Server::builder()
        .add_service(DisciplineServiceServer::from_arc(server.clone()))
        .add_service(ExplicitServiceServer::from_arc(server))
        .serve(address)
        .await?;

    Ok(())
}
