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

//! Requested quantities and execution dependencies.

use crate::aero::SurfaceLoads;
use crate::core::{MatX, Scalar, VecX};
use crate::gradients::GeometryVariable;
use std::collections::BTreeMap;

/// Full geometry remains the default. Half-wing mode never silently falls back.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SolverMode {
    #[default]
    Full,
    HalfWing,
}

/// Scalar outputs, all using the existing full-aircraft normalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Output {
    Lift,
    Drag,
    Side,
    Roll,
    Pitch,
    Yaw,
    NearFieldInducedDrag,
    ViscousDrag,
    FarFieldLift,
    FarFieldSide,
    FarFieldInducedDrag,
    SpanEfficiency,
}

impl Output {
    pub const ALL: [Self; 12] = [
        Self::Lift,
        Self::Drag,
        Self::Side,
        Self::Roll,
        Self::Pitch,
        Self::Yaw,
        Self::NearFieldInducedDrag,
        Self::ViscousDrag,
        Self::FarFieldLift,
        Self::FarFieldSide,
        Self::FarFieldInducedDrag,
        Self::SpanEfficiency,
    ];

    pub(crate) fn far(self) -> bool {
        matches!(
            self,
            Self::FarFieldLift
                | Self::FarFieldSide
                | Self::FarFieldInducedDrag
                | Self::SpanEfficiency
        )
    }
    pub(crate) fn moment(self) -> bool {
        matches!(self, Self::Roll | Self::Pitch | Self::Yaw)
    }
}

/// Independent input of a requested derivative. Angles and rates use the
/// solver's existing radian/normalized-rate conventions. Section twist uses degrees.
#[derive(Clone, Debug, PartialEq)]
pub enum Input {
    Alpha,
    Beta,
    RollRate,
    PitchRate,
    YawRate,
    Mach,
    DynamicPressure,
    MomentReference(usize),
    ReferenceArea,
    ReferenceChord,
    ReferenceSpan,
    Control(String),
    /// Every coordinate of one input surface, identified by its input index.
    MeshNodes(usize),
    Section {
        surface: usize,
        section: usize,
        variable: GeometryVariable,
    },
}

impl Input {
    pub(crate) fn known_zero(&self, output: Output) -> bool {
        matches!(self, Self::DynamicPressure)
            || matches!(self, Self::ReferenceArea) && output == Output::SpanEfficiency
            || matches!(self, Self::ReferenceChord) && !output.moment()
            || matches!(self, Self::ReferenceSpan) && output != Output::SpanEfficiency
    }

    pub(crate) fn geometry(&self) -> bool {
        matches!(self, Self::MeshNodes(_) | Self::Section { .. })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct GradientRequest {
    pub output: Output,
    pub input: Input,
}

/// Optional full-span arrays. Half-wing mode expands only the requested arrays.
#[derive(Clone, Copy, Debug, Default)]
pub struct Details {
    pub circulation: bool,
    pub induced_velocity: bool,
    /// Complete surface and strip load records, as in the existing force API.
    pub surface_loads: bool,
}

/// An empty request performs no aerodynamic work. Derivatives do not implicitly
/// request their scalar outputs. Duplicate requests are evaluated only once.
#[derive(Clone, Debug, Default)]
pub struct AnalysisRequest {
    pub mode: SolverMode,
    pub outputs: Vec<Output>,
    pub gradients: Vec<GradientRequest>,
    pub details: Details,
}

impl AnalysisRequest {
    #[must_use]
    pub fn all_outputs() -> Self {
        Self {
            outputs: Output::ALL.to_vec(),
            ..Self::default()
        }
    }

    pub(crate) fn plan(&self) -> Plan {
        let mut outputs = self.outputs.clone();
        outputs.extend(
            self.gradients
                .iter()
                .filter(|g| !g.input.known_zero(g.output))
                .map(|g| g.output),
        );
        outputs.sort();
        outputs.dedup();
        let near = self.details.surface_loads || outputs.iter().any(|o| !o.far());
        Plan {
            solve: !outputs.is_empty()
                || self.details.circulation
                || self.details.induced_velocity
                || self.details.surface_loads,
            near,
            moments: self.details.surface_loads || outputs.iter().any(|o| o.moment()),
            velocity: near
                || self.details.induced_velocity
                || self
                    .gradients
                    .iter()
                    .any(|g| matches!(g.input, Input::Control(_))),
            far: outputs.iter().any(|o| o.far()),
            far_drag: outputs
                .iter()
                .any(|o| matches!(o, Output::FarFieldInducedDrag | Output::SpanEfficiency)),
            outputs,
        }
    }
}

pub(crate) struct Plan {
    pub solve: bool,
    pub near: bool,
    pub moments: bool,
    pub velocity: bool,
    pub far: bool,
    pub far_drag: bool,
    pub outputs: Vec<Output>,
}

/// Actual work performed during this call (cache hits report zero work).
#[derive(Clone, Debug, Default)]
pub struct Work {
    pub unknowns: usize,
    /// Stored AIC and velocity-operator bytes, excluding the LU and other caches.
    pub operator_bytes: usize,
    pub influence_evaluations: usize,
    pub factorizations: usize,
    pub primal_solves: usize,
    pub velocity_products: usize,
    pub near_field_panels: usize,
    pub moment_panels: usize,
    pub force_derivative_panels: usize,
    pub trefftz_interactions: usize,
    pub adjoint_solves: usize,
    pub mach_interactions: usize,
    pub alpha_derivative_panels: usize,
    pub beta_derivative_panels: usize,
    pub geometry_interactions: usize,
}

#[derive(Clone, Debug)]
pub enum GradientValue {
    Scalar(Scalar),
    /// Node-major x/y/z cotangents in the input surface's original node order.
    MeshNodes(Vec<crate::core::Vec3>),
}

#[derive(Clone, Debug)]
pub struct GradientResult {
    pub request: GradientRequest,
    pub value: GradientValue,
}

/// Omitted quantities are absent, never zero-filled placeholders.
#[derive(Clone, Debug, Default)]
pub struct AnalysisResult {
    pub outputs: BTreeMap<Output, Scalar>,
    pub gradients: Vec<GradientResult>,
    pub circulation: Option<VecX>,
    pub induced_velocity: Option<MatX>,
    pub surface_loads: Option<Vec<SurfaceLoads>>,
    pub work: Work,
}
