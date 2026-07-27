//! Vortex-lattice generation, influence matrices, and primal solver.

mod influence;
mod lattice;
mod types;
mod vlm;

pub use influence::{InfluenceMatrixBuilder, SymmetryConfig, horseshoe_vortex_velocity};
pub use lattice::VortexLattice;
pub(crate) use lattice::chordwise_spacing;
pub use types::{
    ControlSurfaceInfluence, FlightCondition, HorseshoeVortex, SolutionState, Strip,
    SurfaceMeshInfo,
};
pub use vlm::VortexLatticeSolver;
