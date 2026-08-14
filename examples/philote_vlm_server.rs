//! Serves the VLM discipline over gRPC.
//!
//! ```text
//! cargo run --features philote --example philote_vlm_server
//! ```
//!
//! Configure the surfaces through discipline options before `setup`. Both
//! spellings are accepted, so a client written against the v0.2.0 C++
//! discipline works unchanged:
//!
//! ```json
//! {"num_surfaces": 1, "surface_names": ["Wing"],
//!  "surface_nx": [4], "surface_ny": [8]}
//! ```
//! ```json
//! {"num_surfaces": 1, "surface_names": "Wing",
//!  "surface_nx": "4", "surface_ny": "8"}
//! ```

use std::sync::Arc;

use pegasus_vlm::philote::VlmDiscipline;
use philote_mdo::philote_info::{
    discipline_service_server::DisciplineServiceServer,
    explicit_service_server::ExplicitServiceServer,
};
use philote_mdo::server::ExplicitServer;
use tonic::transport::Server;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let address = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "127.0.0.1:50051".to_string())
        .parse()?;
    let server = Arc::new(ExplicitServer::new(VlmDiscipline::default()).with_verbose(true));

    println!("Pegasus VLM discipline listening on {address}");

    Server::builder()
        .add_service(DisciplineServiceServer::from_arc(server.clone()))
        .add_service(ExplicitServiceServer::from_arc(server))
        .serve(address)
        .await?;

    Ok(())
}
