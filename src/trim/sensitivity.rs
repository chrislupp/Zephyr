//! Analytical derivatives of the complete trimmed explicit mapping.

use super::{AeroQuantity, MeshTrimSolution, linearize_aerodynamics};
use crate::core::MatX;
use crate::error::{PegasusError, Result};
use crate::geometry::MeshInput;
use std::collections::BTreeMap;

/// Output/input derivative blocks, flattened in output-major order.
pub type TrimPartials = BTreeMap<(String, String), Vec<f64>>;

impl MeshTrimSolution {
    /// Differentiate requested functions and all solved trim states.
    ///
    /// Input blocks include balance targets and `external_Fx` through
    /// `external_Mz` when a dimensional balance is present. Node blocks can be
    /// omitted when geometry is held fixed. Every block includes trim response.
    pub fn compute_partials(
        &self,
        mesh: &MeshInput,
        outputs: &[AeroQuantity],
        nodes: bool,
    ) -> Result<TrimPartials> {
        self.partials(
            mesh,
            outputs,
            nodes,
            outputs.len() <= self.problem.variables.len(),
        )
    }

    /// Evaluate the same trimmed derivatives using transpose solves for functions.
    pub fn compute_partials_adjoint(
        &self,
        mesh: &MeshInput,
        outputs: &[AeroQuantity],
        nodes: bool,
    ) -> Result<TrimPartials> {
        self.partials(mesh, outputs, nodes, true)
    }

    fn partials(
        &self,
        mesh: &MeshInput,
        outputs: &[AeroQuantity],
        nodes: bool,
        adjoint: bool,
    ) -> Result<TrimPartials> {
        if mesh != &self.mesh {
            return Err(PegasusError::InvalidTrim(
                "derivative mesh must match the converged trim mesh".into(),
            ));
        }
        let quantities: Vec<_> = outputs
            .iter()
            .copied()
            .chain(self.problem.balances.iter().map(|b| b.quantity))
            .collect();
        let mut aero =
            linearize_aerodynamics(mesh, &self.solver, &self.condition, &quantities, nodes)?;
        let fz = self.problem.state_partials(&aero);
        aero.partials.retain(|name, _| !self.problem.owns(name));
        for balance in &self.problem.balances {
            aero.partials
                .insert(balance.target_input.clone(), MatX::zeros(12, 1));
        }
        if self
            .problem
            .balances
            .iter()
            .any(|b| b.quantity.load_index().is_some())
        {
            for q in &AeroQuantity::ALL[6..] {
                aero.partials
                    .insert(format!("external_{}", q.name()), MatX::zeros(12, 1));
            }
        }
        let mut result = BTreeMap::new();
        let lambdas = if adjoint {
            let rhs = MatX::from_fn(self.problem.variables.len(), outputs.len(), |r, c| {
                fz[(outputs[c] as usize, r)]
            });
            Some(
                self.jacobian
                    .transpose()
                    .lu()
                    .solve(&rhs)
                    .ok_or_else(|| PegasusError::TrimFailed("singular trim adjoint".into()))?,
            )
        } else {
            None
        };
        // The small state solve is shared across every output in each input block.
        for (name, fp) in aero.partials {
            let mut rp = MatX::from_fn(self.problem.balances.len(), fp.ncols(), |r, c| {
                fp[(self.problem.balances[r].quantity as usize, c)] / self.problem.balances[r].scale
            });
            for (row, balance) in self.problem.balances.iter().enumerate() {
                if name == balance.target_input {
                    rp[(row, 0)] -= 1.0 / balance.scale;
                }
                if balance.quantity.load_index().is_some()
                    && name == format!("external_{}", balance.quantity.name())
                {
                    rp[(row, 0)] += 1.0 / balance.scale;
                }
            }
            let zp = self
                .factorization
                .solve(&(-&rp))
                .ok_or_else(|| PegasusError::TrimFailed("singular sensitivity system".into()))?;
            let total = &fp + &fz * &zp;
            for (out_col, output) in outputs.iter().enumerate() {
                result.insert(
                    (output.name().into(), name.clone()),
                    (0..total.ncols())
                        .map(|c| {
                            if let Some(lambda) = &lambdas {
                                fp[(*output as usize, c)]
                                    - lambda.column(out_col).dot(&rp.column(c))
                            } else {
                                total[(*output as usize, c)]
                            }
                        })
                        .collect(),
                );
            }
            for (row, variable) in self.problem.variables.iter().enumerate() {
                result.insert(
                    (variable.name.clone(), name.clone()),
                    (0..zp.ncols()).map(|c| zp[(row, c)]).collect(),
                );
            }
            if let Some(mixer) = &self.problem.mixer {
                for (row, surface) in mixer.surfaces.iter().enumerate() {
                    result.insert(
                        (format!("control_{surface}"), name.clone()),
                        (0..zp.ncols())
                            .map(|c| {
                                mixer
                                    .commands
                                    .iter()
                                    .enumerate()
                                    .map(|(col, command)| {
                                        let state = self
                                            .problem
                                            .variables
                                            .iter()
                                            .position(|v| v.name == *command)
                                            .expect("validated command");
                                        mixer.matrix[row][col] * zp[(state, c)]
                                    })
                                    .sum()
                            })
                            .collect(),
                    );
                }
            }
        }
        Ok(result)
    }
}
