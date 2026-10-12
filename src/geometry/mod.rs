//! Aircraft, lifting-surface, body, airfoil, and direct-mesh definitions.

mod aircraft;
mod airfoil;
mod body;
mod control;
mod mesh;
mod section;
mod surface;

pub use aircraft::Aircraft;
pub use airfoil::Airfoil;
pub use body::{Body, BodySection};
pub use control::MeshControlSurface;
pub use mesh::{MeshInput, ReferenceQuantities, SurfaceMesh};
pub use section::{ControlSurfaceDefinition, DesignVariableDefinition, Section};
pub(crate) use surface::SpanwiseLayout;
pub use surface::Surface;
