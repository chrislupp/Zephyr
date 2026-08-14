//! Trim variables, constraints, and common presets.

use crate::core::Scalar;

/// Meaning of one active trim equation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ConstraintType {
    /// Omit this variable and equation.
    #[default]
    Free,
    /// Hold the variable at the target.
    Fixed,
    /// Match lift coefficient.
    Lift,
    /// Match side-force coefficient.
    Side,
    /// Match rolling-moment coefficient.
    RollMoment,
    /// Match pitching-moment coefficient.
    PitchMoment,
    /// Match yawing-moment coefficient.
    YawMoment,
}

/// A target associated with one operating variable.
#[derive(Clone, Copy, Debug, Default)]
pub struct VariableConstraint {
    /// Constraint kind.
    pub kind: ConstraintType,
    /// Desired variable or coefficient value.
    pub target: Scalar,
}

/// State adjusted by the trim solver.
#[derive(Clone, Debug, Default)]
pub struct OperatingVariables {
    pub alpha: Scalar,
    pub beta: Scalar,
    pub roll_rate: Scalar,
    pub pitch_rate: Scalar,
    pub yaw_rate: Scalar,
    /// Ordered control deflections in radians.
    pub control_deflections: Vec<(String, Scalar)>,
}

impl OperatingVariables {
    /// Total state fields, including controls.
    #[must_use]
    pub fn variable_count(&self) -> usize {
        5 + self.control_deflections.len()
    }
}

/// Equations enabled for a trim calculation.
#[derive(Clone, Debug, Default)]
pub struct TrimConstraints {
    pub alpha: VariableConstraint,
    pub beta: VariableConstraint,
    pub roll_rate: VariableConstraint,
    pub pitch_rate: VariableConstraint,
    pub yaw_rate: VariableConstraint,
    /// Ordered control equations.
    pub controls: Vec<(String, VariableConstraint)>,
}

impl TrimConstraints {
    /// Number of active equations.
    #[must_use]
    pub fn constraint_count(&self) -> usize {
        [
            self.alpha,
            self.beta,
            self.roll_rate,
            self.pitch_rate,
            self.yaw_rate,
        ]
        .iter()
        .filter(|constraint| constraint.kind != ConstraintType::Free)
        .count()
            + self
                .controls
                .iter()
                .filter(|(_, constraint)| constraint.kind != ConstraintType::Free)
                .count()
    }
}

/// Constructors for common trim systems.
pub struct TrimPresets;

impl TrimPresets {
    /// Lift target with zero sideslip and zero rates.
    #[must_use]
    pub fn level_flight(target_lift: Scalar) -> TrimConstraints {
        TrimConstraints {
            alpha: VariableConstraint {
                kind: ConstraintType::Lift,
                target: target_lift,
            },
            beta: fixed(0.0),
            roll_rate: fixed(0.0),
            pitch_rate: fixed(0.0),
            yaw_rate: fixed(0.0),
            controls: Vec::new(),
        }
    }

    /// Level flight plus a zero-pitch-moment control equation.
    #[must_use]
    pub fn level_trimmed_flight(
        target_lift: Scalar,
        elevator_name: impl Into<String>,
    ) -> TrimConstraints {
        let mut result = Self::level_flight(target_lift);
        result.controls.push((
            elevator_name.into(),
            VariableConstraint {
                kind: ConstraintType::PitchMoment,
                target: 0.0,
            },
        ));
        result
    }

    /// Steady turn. Bank is retained in the signature for parity.
    #[must_use]
    pub fn steady_turn(target_lift: Scalar, _bank_angle: Scalar) -> TrimConstraints {
        TrimConstraints {
            alpha: VariableConstraint {
                kind: ConstraintType::Lift,
                target: target_lift,
            },
            beta: fixed(0.0),
            roll_rate: fixed(0.0),
            pitch_rate: VariableConstraint::default(),
            yaw_rate: VariableConstraint::default(),
            controls: Vec::new(),
        }
    }

    /// Lift target at fixed sideslip with zero rates.
    #[must_use]
    pub fn steady_sideslip(target_lift: Scalar, sideslip: Scalar) -> TrimConstraints {
        TrimConstraints {
            beta: fixed(sideslip),
            ..Self::level_flight(target_lift)
        }
    }
}

fn fixed(target: Scalar) -> VariableConstraint {
    VariableConstraint {
        kind: ConstraintType::Fixed,
        target,
    }
}
