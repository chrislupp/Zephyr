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

//! Section parameterization to node mesh, with analytical partials.
//!
//! Ports the v0.2.0 C++ `MeshDiscipline`. The span parameterization is
//! `halfspan` times a normalized `eta`, which makes wingspan a design variable
//! with an exact derivative rather than something the caller must rebuild the
//! mesh to explore.

use std::collections::HashMap;

use async_trait::async_trait;
use philote_mdo::{
    ArrayMap, PartialMap, Result, impl_registry,
    registry::VariableRegistry,
    traits::{Discipline, ExplicitDiscipline},
};

use super::{flat, option_bool, option_usize, scalar, scalar_input, sized_input};

/// Explicit discipline mapping section geometry to a structured node mesh.
///
/// Inputs are per-station: `chord`, `twist` (radians), `xle`, `zle` and `eta`,
/// each of length `ny + 1`, plus the scalar `halfspan`. Outputs are the flat
/// `mesh` and the three reference quantities.
#[derive(Debug)]
pub struct MeshDiscipline {
    registry: VariableRegistry,
    nx: usize,
    ny: usize,
    cosine_chordwise_spacing: bool,
    y_duplicated: bool,
    chord_fractions: Vec<f64>,
}

impl Default for MeshDiscipline {
    fn default() -> Self {
        let mut discipline = Self {
            registry: VariableRegistry::new(false),
            nx: 1,
            ny: 1,
            cosine_chordwise_spacing: false,
            y_duplicated: false,
            chord_fractions: Vec::new(),
        };
        discipline.rebuild_chord_fractions();
        discipline
    }
}

impl MeshDiscipline {
    /// A discipline sized for `nx` chordwise and `ny` spanwise panels.
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self {
        let mut discipline = Self {
            nx,
            ny,
            ..Self::default()
        };
        discipline.rebuild_chord_fractions();
        discipline
    }

    /// Number of spanwise stations, one more than the panel count.
    fn station_count(&self) -> usize {
        self.ny + 1
    }

    /// Flattened mesh length: `(nx + 1) * (ny + 1) * 3`.
    fn mesh_len(&self) -> usize {
        (self.nx + 1) * self.station_count() * 3
    }

    /// Chordwise fractions, uniform or cosine-clustered.
    fn rebuild_chord_fractions(&mut self) {
        let count = self.nx + 1;
        self.chord_fractions = (0..count)
            .map(|index| {
                let ratio = index as f64 / self.nx.max(1) as f64;
                if self.cosine_chordwise_spacing {
                    0.5 * (1.0 - (std::f64::consts::PI * ratio).cos())
                } else {
                    ratio
                }
            })
            .collect();
    }

    /// `2` when the surface is mirrored, `1` otherwise.
    fn symmetry_factor(&self) -> f64 {
        if self.y_duplicated { 2.0 } else { 1.0 }
    }

    /// Per-interval arc length, mean chord, and the two reference sums.
    ///
    /// The arc length is measured in the `y`-`z` plane, so dihedral lengthens a
    /// panel rather than being ignored.
    fn planform(
        &self,
        chord: &[f64],
        zle: &[f64],
        eta: &[f64],
        halfspan: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, f64, f64) {
        let mut arc = vec![0.0; self.ny];
        let mut mean_chord = vec![0.0; self.ny];
        let mut d_eta = vec![0.0; self.ny];
        let mut d_zle = vec![0.0; self.ny];
        let mut half_area = 0.0;
        let mut mac_numerator = 0.0;
        for interval in 0..self.ny {
            d_eta[interval] = eta[interval + 1] - eta[interval];
            d_zle[interval] = zle[interval + 1] - zle[interval];
            let dy = d_eta[interval] * halfspan;
            arc[interval] = dy.hypot(d_zle[interval]);
            mean_chord[interval] = 0.5 * (chord[interval] + chord[interval + 1]);
            half_area += mean_chord[interval] * arc[interval];
            mac_numerator += 0.5
                * (chord[interval] * chord[interval] + chord[interval + 1] * chord[interval + 1])
                * arc[interval];
        }
        (arc, mean_chord, d_eta, d_zle, half_area, mac_numerator)
    }
}

impl Discipline for MeshDiscipline {
    impl_registry!(registry);

    fn name(&self) -> &str {
        "ZephyrMesh"
    }

    fn is_differentiable(&self) -> bool {
        true
    }

    fn provides_gradients(&self) -> bool {
        true
    }

    fn initialize(&mut self) -> Result<()> {
        self.add_option("nx", "int")?;
        self.add_option("ny", "int")?;
        self.add_option("cosine_chordwise_spacing", "bool")?;
        self.add_option("y_duplicated", "bool")
    }

    fn set_options(&mut self, options: &HashMap<String, serde_json::Value>) -> Result<()> {
        if let Some(value) = options.get("nx") {
            self.nx = option_usize(value, "nx")?.max(1);
        }
        if let Some(value) = options.get("ny") {
            self.ny = option_usize(value, "ny")?.max(1);
        }
        if let Some(value) = options.get("cosine_chordwise_spacing") {
            self.cosine_chordwise_spacing = option_bool(value);
        }
        if let Some(value) = options.get("y_duplicated") {
            self.y_duplicated = option_bool(value);
        }
        self.rebuild_chord_fractions();
        Ok(())
    }

    fn setup(&mut self) -> Result<()> {
        let stations = self.station_count();
        self.add_input("chord", &[stations], "m")?;
        self.add_input("twist", &[stations], "rad")?;
        self.add_input("xle", &[stations], "m")?;
        self.add_input("zle", &[stations], "m")?;
        self.add_input("halfspan", &[1], "m")?;
        self.add_input("eta", &[stations], "")?;

        self.add_output("mesh", &[self.mesh_len()], "m")?;
        self.add_output("reference_area", &[1], "m**2")?;
        self.add_output("reference_span", &[1], "m")?;
        self.add_output("reference_chord", &[1], "m")
    }

    fn setup_partials(&mut self) -> Result<()> {
        for input in ["chord", "twist", "xle", "zle", "halfspan", "eta"] {
            self.declare_partials("mesh", input)?;
        }
        // The reference quantities do not depend on `xle` or `twist`: neither
        // moves a station in the y-z plane, and the planform integral is taken
        // there.
        for input in ["chord", "zle", "halfspan", "eta"] {
            self.declare_partials("reference_area", input)?;
            self.declare_partials("reference_chord", input)?;
        }
        self.declare_partials("reference_span", "halfspan")?;
        self.declare_partials("reference_span", "eta")
    }
}

#[async_trait]
impl ExplicitDiscipline for MeshDiscipline {
    async fn compute(&self, inputs: &ArrayMap) -> Result<ArrayMap> {
        let stations = self.station_count();
        let chord = sized_input(inputs, "chord", stations)?;
        let twist = sized_input(inputs, "twist", stations)?;
        let xle = sized_input(inputs, "xle", stations)?;
        let zle = sized_input(inputs, "zle", stations)?;
        let eta = sized_input(inputs, "eta", stations)?;
        let halfspan = scalar_input(inputs, "halfspan")?;

        let mut mesh = vec![0.0; self.mesh_len()];
        for station in 0..stations {
            let (sin_twist, cos_twist) = twist[station].sin_cos();
            let y = eta[station] * halfspan;
            for chordwise in 0..=self.nx {
                let fraction = self.chord_fractions[chordwise];
                let base = (chordwise * stations + station) * 3;
                mesh[base] = xle[station] + fraction * chord[station] * cos_twist;
                mesh[base + 1] = y;
                mesh[base + 2] = zle[station] - fraction * chord[station] * sin_twist;
            }
        }

        let (_, _, _, _, half_area, mac_numerator) = self.planform(&chord, &zle, &eta, halfspan);
        let symmetry = self.symmetry_factor();

        let mut outputs = HashMap::new();
        outputs.insert("mesh".to_string(), flat(mesh));
        outputs.insert("reference_area".to_string(), scalar(symmetry * half_area));
        outputs.insert(
            "reference_span".to_string(),
            scalar(symmetry * (eta[stations - 1] - eta[0]) * halfspan),
        );
        outputs.insert(
            "reference_chord".to_string(),
            scalar(if half_area > 0.0 {
                mac_numerator / half_area
            } else {
                0.0
            }),
        );
        Ok(outputs)
    }

    async fn compute_partials(&self, inputs: &ArrayMap) -> Result<PartialMap> {
        let stations = self.station_count();
        let chord = sized_input(inputs, "chord", stations)?;
        let twist = sized_input(inputs, "twist", stations)?;
        let zle = sized_input(inputs, "zle", stations)?;
        let eta = sized_input(inputs, "eta", stations)?;
        let halfspan = scalar_input(inputs, "halfspan")?;

        let mesh_len = self.mesh_len();
        let jacobian = mesh_len * stations;
        let mut d_chord = vec![0.0; jacobian];
        let mut d_twist = vec![0.0; jacobian];
        let mut d_xle = vec![0.0; jacobian];
        let mut d_zle = vec![0.0; jacobian];
        let mut d_eta = vec![0.0; jacobian];
        let mut d_halfspan = vec![0.0; mesh_len];

        // Each mesh coordinate depends on exactly one station, so every row of
        // these Jacobians has a single nonzero.
        for station in 0..stations {
            let (sin_twist, cos_twist) = twist[station].sin_cos();
            for chordwise in 0..=self.nx {
                let fraction = self.chord_fractions[chordwise];
                let base = (chordwise * stations + station) * 3;

                let row_x = base * stations;
                d_chord[row_x + station] = fraction * cos_twist;
                d_twist[row_x + station] = -fraction * chord[station] * sin_twist;
                d_xle[row_x + station] = 1.0;

                let row_y = (base + 1) * stations;
                d_eta[row_y + station] = halfspan;
                d_halfspan[base + 1] = eta[station];

                let row_z = (base + 2) * stations;
                d_chord[row_z + station] = -fraction * sin_twist;
                d_twist[row_z + station] = -fraction * chord[station] * cos_twist;
                d_zle[row_z + station] = 1.0;
            }
        }

        let (arc, mean_chord, d_eta_seg, d_zle_seg, half_area, mac_numerator) =
            self.planform(&chord, &zle, &eta, halfspan);
        let symmetry = self.symmetry_factor();

        let mut area_d_chord = vec![0.0; stations];
        let mut area_d_zle = vec![0.0; stations];
        let mut area_d_eta = vec![0.0; stations];
        let mut area_d_halfspan = 0.0;

        let mut mac_d_chord = vec![0.0; stations];
        let mut mac_d_zle = vec![0.0; stations];
        let mut mac_d_eta = vec![0.0; stations];
        let mut mac_d_halfspan = 0.0;

        let denominator = half_area * half_area;
        for interval in 0..self.ny {
            area_d_chord[interval] += symmetry * 0.5 * arc[interval];
            area_d_chord[interval + 1] += symmetry * 0.5 * arc[interval];

            if arc[interval] <= 0.0 {
                continue;
            }
            // d(arc)/d(halfspan), d(arc)/d(eta) and d(arc)/d(zle) for this
            // interval; the eta terms carry halfspan squared because eta is
            // normalized.
            let d_arc_d_halfspan =
                d_eta_seg[interval] * d_eta_seg[interval] * halfspan / arc[interval];
            let eta_factor = d_eta_seg[interval] * halfspan * halfspan / arc[interval];
            let d_arc_d_zle_left = -d_zle_seg[interval] / arc[interval];
            let d_arc_d_zle_right = d_zle_seg[interval] / arc[interval];

            area_d_halfspan += symmetry * mean_chord[interval] * d_arc_d_halfspan;
            area_d_eta[interval] += symmetry * mean_chord[interval] * (-eta_factor);
            area_d_eta[interval + 1] += symmetry * mean_chord[interval] * eta_factor;
            area_d_zle[interval] += symmetry * mean_chord[interval] * d_arc_d_zle_left;
            area_d_zle[interval + 1] += symmetry * mean_chord[interval] * d_arc_d_zle_right;

            if half_area <= 0.0 {
                continue;
            }
            // Quotient rule on MAC = numerator / half_area. The symmetry factor
            // cancels: it multiplies numerator and denominator alike.
            let numerator_coefficient = 0.5
                * (chord[interval] * chord[interval] + chord[interval + 1] * chord[interval + 1]);
            let denominator_coefficient = mean_chord[interval];
            let accumulate = |slot: &mut f64, d_arc: f64| {
                *slot += (numerator_coefficient * d_arc * half_area
                    - mac_numerator * denominator_coefficient * d_arc)
                    / denominator;
            };
            accumulate(&mut mac_d_halfspan, d_arc_d_halfspan);
            accumulate(&mut mac_d_eta[interval], -eta_factor);
            accumulate(&mut mac_d_eta[interval + 1], eta_factor);
            accumulate(&mut mac_d_zle[interval], d_arc_d_zle_left);
            accumulate(&mut mac_d_zle[interval + 1], d_arc_d_zle_right);
        }

        if half_area > 0.0 {
            for interval in 0..self.ny {
                let d_numerator_left = chord[interval] * arc[interval];
                let d_numerator_right = chord[interval + 1] * arc[interval];
                let d_denominator = 0.5 * arc[interval];
                mac_d_chord[interval] +=
                    (d_numerator_left * half_area - mac_numerator * d_denominator) / denominator;
                mac_d_chord[interval + 1] +=
                    (d_numerator_right * half_area - mac_numerator * d_denominator) / denominator;
            }
        }

        let mut span_d_eta = vec![0.0; stations];
        span_d_eta[0] = -symmetry * halfspan;
        span_d_eta[stations - 1] += symmetry * halfspan;

        let mut partials: PartialMap = HashMap::new();
        let mut insert = |output: &str, input: &str, values: Vec<f64>| {
            partials.insert((output.to_string(), input.to_string()), flat(values));
        };
        insert("mesh", "chord", d_chord);
        insert("mesh", "twist", d_twist);
        insert("mesh", "xle", d_xle);
        insert("mesh", "zle", d_zle);
        insert("mesh", "eta", d_eta);
        insert("mesh", "halfspan", d_halfspan);

        insert("reference_area", "chord", area_d_chord);
        insert("reference_area", "zle", area_d_zle);
        insert("reference_area", "eta", area_d_eta);
        insert("reference_area", "halfspan", vec![area_d_halfspan]);

        insert("reference_chord", "chord", mac_d_chord);
        insert("reference_chord", "zle", mac_d_zle);
        insert("reference_chord", "eta", mac_d_eta);
        insert("reference_chord", "halfspan", vec![mac_d_halfspan]);

        insert(
            "reference_span",
            "halfspan",
            vec![symmetry * (eta[stations - 1] - eta[0])],
        );
        insert("reference_span", "eta", span_d_eta);
        Ok(partials)
    }
}
