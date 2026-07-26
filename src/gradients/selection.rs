//! Explicit output selection for avoiding unnecessary adjoint work.

/// Scalar aerodynamic objective.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum AeroOutput {
    Lift,
    Drag,
    Side,
    Roll,
    Pitch,
    Yaw,
    InducedDrag,
}

/// Flags controlling the gradient outputs to compute.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OutputSelection {
    pub lift: bool,
    pub drag: bool,
    pub side: bool,
    pub roll: bool,
    pub pitch: bool,
    pub yaw: bool,
    pub induced_drag: bool,
}

impl Default for OutputSelection {
    fn default() -> Self {
        Self::all()
    }
}

impl OutputSelection {
    #[must_use]
    pub const fn all() -> Self {
        Self {
            lift: true,
            drag: true,
            side: true,
            roll: true,
            pitch: true,
            yaw: true,
            induced_drag: true,
        }
    }

    #[must_use]
    pub const fn none() -> Self {
        Self {
            lift: false,
            drag: false,
            side: false,
            roll: false,
            pitch: false,
            yaw: false,
            induced_drag: false,
        }
    }

    #[must_use]
    pub const fn only(output: AeroOutput) -> Self {
        let mut selection = Self::none();
        match output {
            AeroOutput::Lift => selection.lift = true,
            AeroOutput::Drag => selection.drag = true,
            AeroOutput::Side => selection.side = true,
            AeroOutput::Roll => selection.roll = true,
            AeroOutput::Pitch => selection.pitch = true,
            AeroOutput::Yaw => selection.yaw = true,
            AeroOutput::InducedDrag => selection.induced_drag = true,
        }
        selection
    }

    #[must_use]
    pub const fn contains(self, output: AeroOutput) -> bool {
        match output {
            AeroOutput::Lift => self.lift,
            AeroOutput::Drag => self.drag,
            AeroOutput::Side => self.side,
            AeroOutput::Roll => self.roll,
            AeroOutput::Pitch => self.pitch,
            AeroOutput::Yaw => self.yaw,
            AeroOutput::InducedDrag => self.induced_drag,
        }
    }
}
