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

use pegasus_vlm::philote::MeshDiscipline;
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
        .unwrap_or_else(|| "127.0.0.1:50052".to_string())
        .parse()?;
    let server = Arc::new(ExplicitServer::new(MeshDiscipline::default()).with_verbose(true));

    println!("Pegasus mesh discipline listening on {address}");

    Server::builder()
        .add_service(DisciplineServiceServer::from_arc(server.clone()))
        .add_service(ExplicitServiceServer::from_arc(server))
        .serve(address)
        .await?;

    Ok(())
}
