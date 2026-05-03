/*
    Project Pegasus
    Version 0.0.0

    Copyright 2024, Christopher A. Lupp

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "pegasus/gradients/mesh_sensitivity.h"

#include "pegasus/geometry/aircraft.h"
#include "pegasus/solver/vortex_lattice.h"

#include <cmath>

namespace pegasus::gradients {

// ============================================================================
// Mapping Construction
// ============================================================================

SectionStripMapping BuildSectionStripMapping(const Aircraft& aircraft,
                                             const VortexLattice& mesh,
                                             Index surface_idx) {
    SectionStripMapping mapping;

    // Handle y-duplicated surfaces
    const auto& surface_info = mesh.GetSurfaceInfo(surface_idx);
    Index original_surface_idx = surface_info.is_y_duplicated
                                     ? surface_info.original_surface_index
                                     : surface_idx;

    // Get the original surface from aircraft
    // Note: For y-duplicated meshes, we need to find the original aircraft
    // surface
    const Surface* surface = nullptr;
    Index aircraft_surface_count = 0;
    for (Index i = 0; i < aircraft.GetNumSurfaces(); ++i) {
        if (aircraft_surface_count == original_surface_idx) {
            surface = &aircraft.GetSurface(i);
            break;
        }
        aircraft_surface_count++;
        if (aircraft.GetSurface(i).IsYDuplicated()) {
            aircraft_surface_count++;  // Y-duplicate adds another mesh surface
        }
    }

    if (!surface) {
        return mapping;  // Return empty mapping if surface not found
    }

    mapping.surface_index = surface_idx;
    mapping.num_sections = surface->GetNumSections();
    mapping.first_strip = surface_info.first_strip_index;
    mapping.num_strips = surface_info.num_strips;

    // Build the strip indices at each section
    // This replicates OptVL's IPTLOC computation
    mapping.strip_at_section.resize(mapping.num_sections);

    // Compute section arc lengths (y-z distance)
    std::vector<double> arc_lengths(mapping.num_sections);
    arc_lengths[0] = 0.0;
    for (Index isec = 1; isec < mapping.num_sections; ++isec) {
        const Section& prev = surface->GetSection(isec - 1);
        const Section& curr = surface->GetSection(isec);
        Vec3 le_prev = prev.GetLeadingEdgePosition();
        Vec3 le_curr = curr.GetLeadingEdgePosition();
        double dy = le_curr.y() - le_prev.y();
        double dz = le_curr.z() - le_prev.z();
        arc_lengths[isec] =
            arc_lengths[isec - 1] + std::sqrt(dy * dy + dz * dz);
    }

    // Compute number of spans per section interval
    std::vector<Index> spans_per_section(mapping.num_sections - 1);
    Index total_spans = 0;
    for (Index isec = 0; isec < mapping.num_sections - 1; ++isec) {
        Index section_panels =
            surface->GetSection(isec).GetSpanwisePanelCount();
        if (section_panels == 0) {
            section_panels = surface->GetSpanwisePanelCount();
            if (section_panels == 0) {
                section_panels = 10;
            }
            // Distribute proportionally by arc length
            double arc_fraction = (arc_lengths[isec + 1] - arc_lengths[isec]) /
                                  arc_lengths.back();
            section_panels = std::max(
                Index(1), static_cast<Index>(section_panels * arc_fraction));
        }
        spans_per_section[isec] = section_panels;
        total_spans += section_panels;
    }

    // Build strip indices at each section boundary
    Index cumulative_strips = 0;
    mapping.strip_at_section[0] = 0;
    for (Index isec = 0; isec < mapping.num_sections - 1; ++isec) {
        cumulative_strips += spans_per_section[isec];
        mapping.strip_at_section[isec + 1] = cumulative_strips;
    }

    return mapping;
}

MeshArrays ExtractMeshArrays(const VortexLattice& mesh) {
    MeshArrays arrays;

    Index num_vortices = mesh.GetNumVortices();
    arrays.num_vortices = num_vortices;

    arrays.rv1.resize(3 * num_vortices);
    arrays.rv2.resize(3 * num_vortices);
    arrays.rc.resize(3 * num_vortices);
    arrays.normal.resize(3 * num_vortices);
    arrays.chord.resize(num_vortices);
    arrays.surface_idx.resize(num_vortices);
    arrays.component_idx.resize(num_vortices);

    for (Index i = 0; i < num_vortices; ++i) {
        const HorseshoeVortex& vor = mesh.GetVortex(i);

        arrays.rv1[3 * i + 0] = vor.bound_vortex_left.x();
        arrays.rv1[3 * i + 1] = vor.bound_vortex_left.y();
        arrays.rv1[3 * i + 2] = vor.bound_vortex_left.z();

        arrays.rv2[3 * i + 0] = vor.bound_vortex_right.x();
        arrays.rv2[3 * i + 1] = vor.bound_vortex_right.y();
        arrays.rv2[3 * i + 2] = vor.bound_vortex_right.z();

        arrays.rc[3 * i + 0] = vor.control_point.x();
        arrays.rc[3 * i + 1] = vor.control_point.y();
        arrays.rc[3 * i + 2] = vor.control_point.z();

        arrays.normal[3 * i + 0] = vor.normal_at_control_point.x();
        arrays.normal[3 * i + 1] = vor.normal_at_control_point.y();
        arrays.normal[3 * i + 2] = vor.normal_at_control_point.z();

        arrays.chord[i] = vor.local_chord;
        arrays.surface_idx[i] = static_cast<int>(vor.surface_index);
        arrays.component_idx[i] = static_cast<int>(vor.component_index);
    }

    return arrays;
}

void BuildDefaultChordwiseFractions(Index num_chordwise,
                                    std::vector<double>& xvr_out,
                                    std::vector<double>& xcp_out,
                                    SpacingType spacing) {
    xvr_out.resize(num_chordwise);
    xcp_out.resize(num_chordwise);

    // Use the same spacing as the mesh generator (GenerateCosineSpacing for
    // panel edges)
    std::vector<Scalar> x_pt =
        GenerateCosineSpacing(num_chordwise + 1, spacing);

    for (Index i = 0; i < num_chordwise; ++i) {
        double x1 = x_pt[i];
        double x2 = x_pt[i + 1];
        double dx = x2 - x1;

        // Standard VLM: vortex at 1/4 panel chord, control at 3/4
        xvr_out[i] = x1 + 0.25 * dx;
        xcp_out[i] = x1 + 0.75 * dx;
    }
}

// ============================================================================
// Section Gradient Computation
// ============================================================================

SectionGradients ComputeSectionGradientFromMesh(
    const MeshGradients& mesh_grads, const SectionStripMapping& mapping,
    Index section_idx, const VortexLattice& mesh,
    const std::vector<double>& xvr_fractions,
    const std::vector<double>& xcp_fractions) {
    SectionGradients result;

    // Get the range of strips affected by this section interval
    // Section section_idx affects the interval [section_idx, section_idx+1]
    auto [first_strip, num_strips] = mapping.GetStripRange(section_idx);

    if (num_strips == 0) {
        return result;
    }

    // Get surface info
    const auto& surface_info = mesh.GetSurfaceInfo(mapping.surface_index);
    Index num_chordwise = surface_info.num_chordwise;

    // Accumulators for left and right section contributions
    double d_le_left[3] = {0.0, 0.0, 0.0};
    double d_le_right[3] = {0.0, 0.0, 0.0};
    double d_chord_left = 0.0;
    double d_chord_right = 0.0;
    double d_twist_left = 0.0;
    double d_twist_right = 0.0;

    // Loop over strips in this section interval
    for (Index istrip = 0; istrip < num_strips; ++istrip) {
        Index strip_idx = first_strip + istrip;
        const Strip& strip = mesh.GetStrip(strip_idx);

        // Compute interpolation fractions
        double f1, f2;
        ComputeStripEdgeFractions(istrip, num_strips, &f1, &f2);
        double fc = ComputeStripInterpolationFraction(istrip, num_strips);

        // Get Trefftz plane normal for twist gradient computation
        double trefftz_ny = strip.trefftz_normal_y;
        double trefftz_nz = strip.trefftz_normal_z;

        // Loop over vortex elements in this strip
        for (Index ichord = 0; ichord < num_chordwise; ++ichord) {
            Index vor_idx = strip.first_vortex_index + ichord;
            const HorseshoeVortex& vor = mesh.GetVortex(vor_idx);

            // Get chordwise fractions
            double xvr = (ichord < static_cast<Index>(xvr_fractions.size()))
                             ? xvr_fractions[ichord]
                             : 0.25;
            double xcp = (ichord < static_cast<Index>(xcp_fractions.size()))
                             ? xcp_fractions[ichord]
                             : 0.75;

            // Get mesh gradients for this vortex
            const double* d_rv1 = &mesh_grads.d_rv1[3 * vor_idx];
            const double* d_rv2 = &mesh_grads.d_rv2[3 * vor_idx];
            const double* d_rc = &mesh_grads.d_rc[3 * vor_idx];
            const double* d_normal = &mesh_grads.d_normal[3 * vor_idx];

            // -----------------------------------------------------------------
            // Leading edge gradient contribution
            // -----------------------------------------------------------------
            // RV1 contribution (left edge uses RLE1)
            double d_rle1[3] = {d_rv1[0], d_rv1[1], d_rv1[2]};
            PropagateLeadingEdgeGradient(d_rle1, f1, d_le_left, d_le_right);

            // RV2 contribution (right edge uses RLE2)
            double d_rle2[3] = {d_rv2[0], d_rv2[1], d_rv2[2]};
            PropagateLeadingEdgeGradient(d_rle2, f2, d_le_left, d_le_right);

            // RC contribution (center uses RLE)
            double d_rle[3] = {d_rc[0], d_rc[1], d_rc[2]};
            PropagateLeadingEdgeGradient(d_rle, fc, d_le_left, d_le_right);

            // -----------------------------------------------------------------
            // Chord gradient contribution
            // -----------------------------------------------------------------
            PropagateChordGradient(d_rv1, d_rv2, d_rc, xvr, xcp, f1, f2, fc,
                                   &d_chord_left, &d_chord_right);

            // -----------------------------------------------------------------
            // Twist gradient contribution
            // -----------------------------------------------------------------
            // Compute bound direction for this vortex
            double bound_dir[3];
            double dx = vor.bound_vortex_right.x() - vor.bound_vortex_left.x();
            double dy = vor.bound_vortex_right.y() - vor.bound_vortex_left.y();
            double dz = vor.bound_vortex_right.z() - vor.bound_vortex_left.z();
            double bound_len = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (bound_len > 1e-12) {
                bound_dir[0] = dx / bound_len;
                bound_dir[1] = dy / bound_len;
                bound_dir[2] = dz / bound_len;
            } else {
                bound_dir[0] = 0.0;
                bound_dir[1] = 1.0;
                bound_dir[2] = 0.0;
            }

            // Get camber slope at control point
            double camber_slope = vor.camber_slope_at_control_point;

            // Compute twist gradient contribution
            double d_twist_contribution =
                PropagateTwistGradient(d_normal, strip.twist_angle, trefftz_ny,
                                       trefftz_nz, bound_dir, camber_slope);

            // Distribute to left and right sections based on interpolation
            // The twist at strip center is interpolated from sections:
            // twist(strip) = atan2(CHSIN, CHCOS) where
            //   CHSIN = (1-fc)*CHORDL*sin(AINCL) + fc*CHORDR*sin(AINCR)
            //   CHCOS = (1-fc)*CHORDL*cos(AINCL) + fc*CHORDR*cos(AINCR)
            // For small angles, this is approximately linear:
            //   twist ≈ (1-fc)*twist_left + fc*twist_right
            d_twist_left += (1.0 - fc) * d_twist_contribution;
            d_twist_right += fc * d_twist_contribution;
        }
    }

    // The gradient computed is for section_idx (left) and section_idx+1 (right)
    // We return gradients for section_idx only; the caller must also query
    // section_idx-1 to get the "right" contribution from the previous interval
    result.d_le_x = d_le_left[0];
    result.d_le_y = d_le_left[1];
    result.d_le_z = d_le_left[2];
    result.d_chord = d_chord_left;
    result.d_twist = d_twist_left;

    // Note: For interior sections, we need to add contributions from both
    // the interval before (where this section is "right") and the interval
    // after (where this section is "left"). The caller should handle this.

    return result;
}

}  // namespace pegasus::gradients
