//! Cached explicit aerodynamic evaluation, optionally containing internal trim.

use crate::geometry::MeshInput;
use crate::solver::{FlightCondition, VortexLatticeSolver};
use crate::trim::MeshTrimSolution;
use philote_mdo::{ArrayMap, PartialMap};

#[derive(Debug)]
pub(super) struct SolvedState {
    pub inputs: ArrayMap,
    pub mesh: MeshInput,
    pub condition: FlightCondition,
    pub solver: VortexLatticeSolver,
    pub trim: Option<MeshTrimSolution>,
    pub partials: Option<PartialMap>,
}

impl SolvedState {
    pub fn solver(&self) -> &VortexLatticeSolver {
        self.trim.as_ref().map_or(&self.solver, |t| &t.solver)
    }
}
