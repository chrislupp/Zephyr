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
