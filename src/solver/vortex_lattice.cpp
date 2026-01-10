// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/solver/vortex_lattice.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "vlm/utils/math_utilities.h"

namespace vlm {

void VortexLattice::Clear() {
  vortices_.clear();
  strips_.clear();
  surface_info_.clear();
  control_influences_.clear();
  control_surface_names_.clear();
}

void VortexLattice::GenerateFromAircraft(const Aircraft& aircraft) {
  Clear();

  // Collect control surface names from the aircraft
  control_surface_names_ = aircraft.GetControlSurfaceNames();

  // Generate mesh for each surface
  for (Index i = 0; i < aircraft.GetNumSurfaces(); ++i) {
    const Surface& surface = aircraft.GetSurface(i);
    GenerateSurfaceMesh(surface, i, aircraft);

    // Handle y-duplication
    if (surface.IsYDuplicated()) {
      GenerateDuplicateSurfaceMesh(
          static_cast<Index>(surface_info_.size()) - 1,
          surface.GetYDuplicate().value());
    }
  }

  // Calculate normal vectors for all elements
  CalculateNormals();
}

void VortexLattice::GenerateSurfaceMesh(const Surface& surface,
                                        Index surface_index,
                                        const Aircraft& aircraft) {
  if (surface.GetNumSections() < 2) {
    throw std::runtime_error(
        "VortexLattice::GenerateSurfaceMesh: Need at least 2 sections per "
        "surface");
  }

  SurfaceMeshInfo info;
  info.name = surface.GetName();
  info.first_vortex_index = static_cast<Index>(vortices_.size());
  info.first_strip_index = static_cast<Index>(strips_.size());
  info.num_chordwise = surface.GetChordwisePanelCount();
  info.is_y_duplicated = false;
  info.original_surface_index = surface_index;

  // Get surface transformation parameters
  Vec3 scale = surface.GetScale();
  Vec3 translation = surface.GetTranslation();
  Scalar angle_offset = surface.GetAngleOffset() * kDegreesToRadians;

  // Calculate arc lengths along the wing trace (y-z plane)
  std::vector<Scalar> arc_lengths(surface.GetNumSections());
  arc_lengths[0] = 0.0;
  for (Index isec = 1; isec < surface.GetNumSections(); ++isec) {
    const Section& prev = surface.GetSection(isec - 1);
    const Section& curr = surface.GetSection(isec);
    Vec3 le_prev = prev.GetLeadingEdgePosition();
    Vec3 le_curr = curr.GetLeadingEdgePosition();
    Scalar dy = le_curr.y() - le_prev.y();
    Scalar dz = le_curr.z() - le_prev.z();
    arc_lengths[isec] = arc_lengths[isec - 1] + std::sqrt(dy * dy + dz * dz);
  }

  // Calculate spanwise panel distribution
  // First determine total spanwise panels based on section intervals
  Index total_spanwise = 0;
  std::vector<Index> span_panels_per_section(surface.GetNumSections() - 1);

  for (Index isec = 0; isec < surface.GetNumSections() - 1; ++isec) {
    const Section& sec = surface.GetSection(isec);
    Index section_panels = sec.GetSpanwisePanelCount();
    if (section_panels == 0) {
      // Use surface default if section doesn't specify
      section_panels = surface.GetSpanwisePanelCount();
      if (section_panels == 0) {
        section_panels = 10;  // Reasonable default
      }
      // Distribute proportionally by arc length
      Scalar arc_fraction =
          (arc_lengths[isec + 1] - arc_lengths[isec]) / arc_lengths.back();
      section_panels =
          std::max(Index(1), static_cast<Index>(section_panels * arc_fraction));
    }
    span_panels_per_section[isec] = section_panels;
    total_spanwise += section_panels;
  }

  info.num_spanwise = total_spanwise;

  // Generate spanwise spacing fractions for each section interval
  std::vector<Scalar> y_panel_edges;  // Panel edge positions
  std::vector<Scalar> y_control;      // Control point positions
  y_panel_edges.push_back(arc_lengths[0]);

  for (Index isec = 0; isec < surface.GetNumSections() - 1; ++isec) {
    const Section& sec = surface.GetSection(isec);
    Index num_span = span_panels_per_section[isec];
    Scalar arc_start = arc_lengths[isec];
    Scalar arc_end = arc_lengths[isec + 1];
    Scalar arc_length = arc_end - arc_start;

    // Get spacing type for this interval
    SpacingType spacing = sec.GetSpanwiseSpacing();
    if (spacing == SpacingType::kNone) {
      spacing = surface.GetSpanwiseSpacing();
    }

    // Generate spacing fractions
    std::vector<Scalar> spacing_fracs =
        GenerateCosineSpacing(2 * num_span + 1, spacing);

    for (Index n = 0; n < num_span; ++n) {
      // Control point at odd indices (2n+1)
      y_control.push_back(arc_start + arc_length * spacing_fracs[2 * n + 1]);
      // Panel edge at even indices (2n+2)
      y_panel_edges.push_back(arc_start +
                              arc_length * spacing_fracs[2 * n + 2]);
    }
  }

  // Generate chordwise spacing
  std::vector<Scalar> x_vortex(info.num_chordwise);
  std::vector<Scalar> x_control(info.num_chordwise);
  CalculateChordwiseSpacing(info.num_chordwise, surface.GetChordwiseSpacing(),
                            1.0, x_vortex, x_control);

  // Track cumulative wetted area
  Scalar wetted_area = 0.0;
  Scalar total_width = 0.0;

  // Generate strips and vortices for each section interval
  Index strip_idx = 0;
  for (Index isec = 0; isec < surface.GetNumSections() - 1; ++isec) {
    const Section& sec_left = surface.GetSection(isec);
    const Section& sec_right = surface.GetSection(isec + 1);

    // Apply scale and translation to section geometry
    Vec3 le_left;
    le_left.x() = scale.x() * sec_left.GetLeadingEdgePosition().x() +
                  translation.x();
    le_left.y() = scale.y() * sec_left.GetLeadingEdgePosition().y() +
                  translation.y();
    le_left.z() = scale.z() * sec_left.GetLeadingEdgePosition().z() +
                  translation.z();

    Vec3 le_right;
    le_right.x() = scale.x() * sec_right.GetLeadingEdgePosition().x() +
                   translation.x();
    le_right.y() = scale.y() * sec_right.GetLeadingEdgePosition().y() +
                   translation.y();
    le_right.z() = scale.z() * sec_right.GetLeadingEdgePosition().z() +
                   translation.z();

    Scalar chord_left = scale.x() * sec_left.GetChord();
    Scalar chord_right = scale.x() * sec_right.GetChord();

    // Section incidence angles (include angle offset)
    Scalar twist_left =
        sec_left.GetTwistAngle() * kDegreesToRadians + angle_offset;
    Scalar twist_right =
        sec_right.GetTwistAngle() * kDegreesToRadians + angle_offset;

    // CL-alpha factors
    Scalar claf_left = sec_left.GetClAlphaFactor();
    Scalar claf_right = sec_right.GetClAlphaFactor();

    // Section width
    Scalar section_width =
        std::sqrt(std::pow(le_right.y() - le_left.y(), 2) +
                  std::pow(le_right.z() - le_left.z(), 2));

    // Leading/trailing edge sweep slopes
    Scalar sweep_le = (section_width > 0)
                          ? (le_right.x() - le_left.x()) / section_width
                          : 0.0;
    Scalar sweep_te = (section_width > 0)
                          ? (le_right.x() + chord_right - le_left.x() -
                             chord_left) /
                                section_width
                          : 0.0;

    // Get control surfaces for this section interval
    const auto& controls_left = sec_left.GetControlSurfaces();
    const auto& controls_right = sec_right.GetControlSurfaces();

    // Precompute sin/cos of twist for interpolation
    Scalar sin_twist_left = chord_left * std::sin(twist_left);
    Scalar cos_twist_left = chord_left * std::cos(twist_left);
    Scalar sin_twist_right = chord_right * std::sin(twist_right);
    Scalar cos_twist_right = chord_right * std::cos(twist_right);

    // Generate strips for this section interval
    Index num_span = span_panels_per_section[isec];
    Scalar arc_start = arc_lengths[isec];
    Scalar arc_end = arc_lengths[isec + 1];

    for (Index ispan = 0; ispan < num_span; ++ispan) {
      Strip strip;

      // Get panel edge and center arc positions
      Scalar y_pt1 = y_panel_edges[strip_idx];
      Scalar y_pt2 = y_panel_edges[strip_idx + 1];
      Scalar y_cp = y_control[strip_idx];

      // Interpolation fractions within this section interval
      Scalar f1 = (y_pt1 - arc_start) / (arc_end - arc_start);
      Scalar f2 = (y_pt2 - arc_start) / (arc_end - arc_start);
      Scalar fc = (y_cp - arc_start) / (arc_end - arc_start);

      // Interpolate leading edge positions
      strip.leading_edge_left = (1.0 - f1) * le_left + f1 * le_right;
      strip.leading_edge_right = (1.0 - f2) * le_left + f2 * le_right;
      strip.leading_edge_center = (1.0 - fc) * le_left + fc * le_right;

      // Interpolate chords
      strip.chord_left = (1.0 - f1) * chord_left + f1 * chord_right;
      strip.chord_right = (1.0 - f2) * chord_left + f2 * chord_right;
      strip.chord_center = (1.0 - fc) * chord_left + fc * chord_right;

      // Strip width
      strip.width = std::abs(f2 - f1) * section_width;

      // Sweep slopes
      strip.leading_edge_sweep_slope = sweep_le;
      strip.trailing_edge_sweep_slope = sweep_te;

      // Interpolate twist angles using chord-weighted sin/cos
      Scalar sin_interp = (1.0 - fc) * sin_twist_left + fc * sin_twist_right;
      Scalar cos_interp = (1.0 - fc) * cos_twist_left + fc * cos_twist_right;
      strip.twist_angle = std::atan2(sin_interp, cos_interp);

      // Left/right edge twist
      Scalar sin_l = (1.0 - f1) * sin_twist_left + f1 * sin_twist_right;
      Scalar cos_l = (1.0 - f1) * cos_twist_left + f1 * cos_twist_right;
      strip.twist_angle_left = std::atan2(sin_l, cos_l);

      Scalar sin_r = (1.0 - f2) * sin_twist_left + f2 * sin_twist_right;
      Scalar cos_r = (1.0 - f2) * cos_twist_left + f2 * cos_twist_right;
      strip.twist_angle_right = std::atan2(sin_r, cos_r);

      // Trefftz plane normal (y-z projection)
      Scalar dy = strip.leading_edge_right.y() - strip.leading_edge_left.y();
      Scalar dz = strip.leading_edge_right.z() - strip.leading_edge_left.z();
      Scalar yz_mag = std::sqrt(dy * dy + dz * dz);
      if (yz_mag > kEpsilon) {
        strip.trefftz_normal_y = -dz / yz_mag;
        strip.trefftz_normal_z = dy / yz_mag;
      }

      // Reference point
      strip.reference_point = strip.leading_edge_center;
      strip.reference_point.x() +=
          spanwise_axis_fraction_ * strip.chord_center;

      // Spanwise unit vector
      if (yz_mag > kEpsilon) {
        strip.spanwise_unit_vector.x() = 0.0;
        strip.spanwise_unit_vector.y() = dy / yz_mag;
        strip.spanwise_unit_vector.z() = dz / yz_mag;
      }

      // Indices
      strip.first_vortex_index = static_cast<Index>(vortices_.size());
      strip.num_vortices = info.num_chordwise;
      strip.surface_index = static_cast<Index>(surface_info_.size());

      // Interpolate viscous polar
      const auto& polar_left = sec_left.GetViscousPolar();
      const auto& polar_right = sec_right.GetViscousPolar();
      for (int l = 0; l < 6; ++l) {
        strip.viscous_polar[l] =
            (1.0 - fc) * polar_left[l] + fc * polar_right[l];
      }
      strip.has_viscous_data = (strip.viscous_polar[3] != 0.0);

      // Accumulate wetted area
      wetted_area += strip.width * strip.chord_center;
      total_width += strip.width;

      strips_.push_back(strip);

      // Generate horseshoe vortex elements for this strip
      Scalar claf_c = (1.0 - fc) * (chord_left / strip.chord_center) *
                          claf_left +
                      fc * (chord_right / strip.chord_center) * claf_right;

      // Recalculate chordwise spacing with local CL-alpha factor
      std::vector<Scalar> x_vr(info.num_chordwise);
      std::vector<Scalar> x_cp(info.num_chordwise);
      CalculateChordwiseSpacing(info.num_chordwise,
                                surface.GetChordwiseSpacing(), claf_c, x_vr,
                                x_cp);

      // Also need panel edge positions for element widths
      std::vector<Scalar> x_pt =
          GenerateCosineSpacing(info.num_chordwise + 1,
                                surface.GetChordwiseSpacing());

      for (Index ivc = 0; ivc < info.num_chordwise; ++ivc) {
        HorseshoeVortex vortex;

        // Bound vortex endpoints (at panel edges, spanwise)
        vortex.bound_vortex_left.x() =
            strip.leading_edge_left.x() + x_vr[ivc] * strip.chord_left;
        vortex.bound_vortex_left.y() = strip.leading_edge_left.y();
        vortex.bound_vortex_left.z() = strip.leading_edge_left.z();

        vortex.bound_vortex_right.x() =
            strip.leading_edge_right.x() + x_vr[ivc] * strip.chord_right;
        vortex.bound_vortex_right.y() = strip.leading_edge_right.y();
        vortex.bound_vortex_right.z() = strip.leading_edge_right.z();

        // Bound vortex center
        vortex.bound_vortex_center.x() =
            strip.leading_edge_center.x() + x_vr[ivc] * strip.chord_center;
        vortex.bound_vortex_center.y() = strip.leading_edge_center.y();
        vortex.bound_vortex_center.z() = strip.leading_edge_center.z();

        // Control point
        vortex.control_point.x() =
            strip.leading_edge_center.x() + x_cp[ivc] * strip.chord_center;
        vortex.control_point.y() = strip.leading_edge_center.y();
        vortex.control_point.z() = strip.leading_edge_center.z();

        // Element properties
        vortex.chordwise_width = (x_pt[ivc + 1] - x_pt[ivc]) * strip.chord_center;
        vortex.local_chord = strip.chord_center;
        vortex.strip_index = static_cast<Index>(strips_.size()) - 1;
        vortex.surface_index = static_cast<Index>(surface_info_.size());
        vortex.component_index = surface.GetComponentIndex();

        // Camber slopes from airfoil (if available)
        Scalar slope_cp = 0.0;
        Scalar slope_vr = 0.0;
        if (sec_left.HasAirfoil()) {
          const Airfoil& af = sec_left.GetAirfoil();
          slope_cp = (1.0 - fc) * af.EvaluateCamberSlope(x_cp[ivc]);
          slope_vr = (1.0 - fc) * af.EvaluateCamberSlope(x_vr[ivc]);
        }
        if (sec_right.HasAirfoil()) {
          const Airfoil& af = sec_right.GetAirfoil();
          slope_cp += fc * af.EvaluateCamberSlope(x_cp[ivc]);
          slope_vr += fc * af.EvaluateCamberSlope(x_vr[ivc]);
        }
        vortex.camber_slope_at_control_point = slope_cp;
        vortex.camber_slope_at_vortex = slope_vr;

        // Flags
        vortex.affected_by_freestream = surface.IsAffectedByAlphaBeta();
        vortex.enforce_flow_tangency = surface.HasWake();

        vortices_.push_back(vortex);

        // Control surface influences
        std::vector<ControlSurfaceInfluence> influences;
        for (const auto& ctrl : controls_left) {
          // Find matching control on right section
          const ControlSurfaceDefinition* ctrl_right = nullptr;
          for (const auto& cr : controls_right) {
            if (cr.name == ctrl.name) {
              ctrl_right = &cr;
              break;
            }
          }

          if (ctrl_right != nullptr) {
            ControlSurfaceInfluence influence;
            influence.name = ctrl.name;

            // Interpolate hinge position
            Scalar x_hinge_left = ctrl.hinge_x_fraction * chord_left;
            Scalar x_hinge_right = ctrl_right->hinge_x_fraction * chord_right;
            Scalar x_hinge = (1.0 - fc) * x_hinge_left + fc * x_hinge_right;

            // Calculate gain (fraction of element on control surface)
            Scalar x_le_ctrl, x_te_ctrl;
            if (x_hinge >= 0) {
              // Trailing edge control surface
              x_le_ctrl = x_hinge / strip.chord_center;
              x_te_ctrl = 1.0;
            } else {
              // Leading edge control surface
              x_le_ctrl = 0.0;
              x_te_ctrl = -x_hinge / strip.chord_center;
            }

            Scalar frac_le = (x_le_ctrl - x_pt[ivc]) / (x_pt[ivc + 1] - x_pt[ivc]);
            Scalar frac_te = (x_te_ctrl - x_pt[ivc]) / (x_pt[ivc + 1] - x_pt[ivc]);
            frac_le = std::max(0.0, std::min(1.0, frac_le));
            frac_te = std::max(0.0, std::min(1.0, frac_te));

            Scalar gain = (1.0 - fc) * ctrl.gain + fc * ctrl_right->gain;
            influence.gain = gain * (frac_te - frac_le);

            // Hinge axis (interpolate and normalize)
            Vec3 hinge_vec = (1.0 - fc) * ctrl.hinge_vector +
                             fc * ctrl_right->hinge_vector;
            Scalar hinge_mag = hinge_vec.norm();
            if (hinge_mag > kEpsilon) {
              influence.hinge_axis = hinge_vec / hinge_mag;
            } else {
              // Default hinge vector along span
              influence.hinge_axis = strip.spanwise_unit_vector;
            }

            // Hinge point
            influence.hinge_point = strip.leading_edge_center;
            influence.hinge_point.x() += std::abs(x_hinge);

            influence.reflection_sign = ctrl.sign;

            if (std::abs(influence.gain) > kEpsilon) {
              influences.push_back(influence);
            }
          }
        }
        control_influences_.push_back(influences);
      }

      ++strip_idx;
    }
  }

  info.num_vortices =
      static_cast<Index>(vortices_.size()) - info.first_vortex_index;
  info.num_strips =
      static_cast<Index>(strips_.size()) - info.first_strip_index;
  info.wetted_area = wetted_area;
  info.average_chord = (total_width > 0) ? wetted_area / total_width : 0.0;

  surface_info_.push_back(info);
}

void VortexLattice::GenerateDuplicateSurfaceMesh(Index original_surface_index,
                                                 Scalar y_duplicate_plane) {
  const SurfaceMeshInfo& orig = surface_info_[original_surface_index];

  SurfaceMeshInfo info;
  info.name = orig.name + " (yduplicate)";
  info.first_vortex_index = static_cast<Index>(vortices_.size());
  info.first_strip_index = static_cast<Index>(strips_.size());
  info.num_chordwise = orig.num_chordwise;
  info.num_spanwise = orig.num_spanwise;
  info.is_y_duplicated = true;
  info.original_surface_index = original_surface_index;
  info.wetted_area = orig.wetted_area;
  info.average_chord = orig.average_chord;

  Scalar y_offset = 2.0 * y_duplicate_plane;

  // Pre-allocate to avoid reallocation during the loop
  // This prevents iterator/reference invalidation when we push_back
  strips_.reserve(strips_.size() + orig.num_strips);
  vortices_.reserve(vortices_.size() + orig.num_vortices);
  control_influences_.reserve(control_influences_.size() + orig.num_vortices);

  // Create mirrored strips (with edges swapped to maintain circulation sign)
  for (Index j = 0; j < orig.num_strips; ++j) {
    const Strip& orig_strip = strips_[orig.first_strip_index + j];
    Strip new_strip;

    // Swap left/right and reflect y
    new_strip.leading_edge_left.x() = orig_strip.leading_edge_right.x();
    new_strip.leading_edge_left.y() =
        -orig_strip.leading_edge_right.y() + y_offset;
    new_strip.leading_edge_left.z() = orig_strip.leading_edge_right.z();

    new_strip.leading_edge_right.x() = orig_strip.leading_edge_left.x();
    new_strip.leading_edge_right.y() =
        -orig_strip.leading_edge_left.y() + y_offset;
    new_strip.leading_edge_right.z() = orig_strip.leading_edge_left.z();

    new_strip.leading_edge_center.x() = orig_strip.leading_edge_center.x();
    new_strip.leading_edge_center.y() =
        -orig_strip.leading_edge_center.y() + y_offset;
    new_strip.leading_edge_center.z() = orig_strip.leading_edge_center.z();

    new_strip.chord_left = orig_strip.chord_right;
    new_strip.chord_right = orig_strip.chord_left;
    new_strip.chord_center = orig_strip.chord_center;

    new_strip.width = orig_strip.width;
    new_strip.leading_edge_sweep_slope = -orig_strip.leading_edge_sweep_slope;
    new_strip.trailing_edge_sweep_slope = -orig_strip.trailing_edge_sweep_slope;

    new_strip.twist_angle = orig_strip.twist_angle;
    new_strip.twist_angle_left = orig_strip.twist_angle_right;
    new_strip.twist_angle_right = orig_strip.twist_angle_left;

    new_strip.trefftz_normal_y = -orig_strip.trefftz_normal_y;
    new_strip.trefftz_normal_z = orig_strip.trefftz_normal_z;

    new_strip.spanwise_unit_vector.x() = orig_strip.spanwise_unit_vector.x();
    new_strip.spanwise_unit_vector.y() = -orig_strip.spanwise_unit_vector.y();
    new_strip.spanwise_unit_vector.z() = orig_strip.spanwise_unit_vector.z();

    new_strip.reference_point.x() = orig_strip.reference_point.x();
    new_strip.reference_point.y() =
        -orig_strip.reference_point.y() + y_offset;
    new_strip.reference_point.z() = orig_strip.reference_point.z();

    new_strip.first_vortex_index = static_cast<Index>(vortices_.size());
    new_strip.num_vortices = orig_strip.num_vortices;
    new_strip.surface_index = static_cast<Index>(surface_info_.size());

    new_strip.viscous_polar = orig_strip.viscous_polar;
    new_strip.has_viscous_data = orig_strip.has_viscous_data;
    new_strip.is_active = orig_strip.is_active;

    strips_.push_back(new_strip);

    // Create mirrored vortex elements
    for (Index k = 0; k < orig_strip.num_vortices; ++k) {
      const HorseshoeVortex& orig_vor =
          vortices_[orig_strip.first_vortex_index + k];
      HorseshoeVortex new_vor;

      // Swap left/right and reflect y
      new_vor.bound_vortex_left.x() = orig_vor.bound_vortex_right.x();
      new_vor.bound_vortex_left.y() =
          -orig_vor.bound_vortex_right.y() + y_offset;
      new_vor.bound_vortex_left.z() = orig_vor.bound_vortex_right.z();

      new_vor.bound_vortex_right.x() = orig_vor.bound_vortex_left.x();
      new_vor.bound_vortex_right.y() =
          -orig_vor.bound_vortex_left.y() + y_offset;
      new_vor.bound_vortex_right.z() = orig_vor.bound_vortex_left.z();

      new_vor.bound_vortex_center.x() = orig_vor.bound_vortex_center.x();
      new_vor.bound_vortex_center.y() =
          -orig_vor.bound_vortex_center.y() + y_offset;
      new_vor.bound_vortex_center.z() = orig_vor.bound_vortex_center.z();

      new_vor.control_point.x() = orig_vor.control_point.x();
      new_vor.control_point.y() = -orig_vor.control_point.y() + y_offset;
      new_vor.control_point.z() = orig_vor.control_point.z();

      new_vor.camber_slope_at_control_point =
          orig_vor.camber_slope_at_control_point;
      new_vor.camber_slope_at_vortex = orig_vor.camber_slope_at_vortex;
      new_vor.chordwise_width = orig_vor.chordwise_width;
      new_vor.local_chord = orig_vor.local_chord;
      new_vor.strip_index = static_cast<Index>(strips_.size()) - 1;
      new_vor.surface_index = static_cast<Index>(surface_info_.size());
      new_vor.component_index = orig_vor.component_index;
      new_vor.affected_by_freestream = orig_vor.affected_by_freestream;
      new_vor.enforce_flow_tangency = orig_vor.enforce_flow_tangency;

      vortices_.push_back(new_vor);

      // Mirror control surface influences
      const auto& orig_influences =
          control_influences_[orig_strip.first_vortex_index + k];
      std::vector<ControlSurfaceInfluence> new_influences;
      for (const auto& orig_inf : orig_influences) {
        ControlSurfaceInfluence new_inf;
        new_inf.name = orig_inf.name;
        new_inf.gain = -orig_inf.gain * orig_inf.reflection_sign;
        new_inf.hinge_axis.x() = orig_inf.hinge_axis.x();
        new_inf.hinge_axis.y() = -orig_inf.hinge_axis.y();
        new_inf.hinge_axis.z() = orig_inf.hinge_axis.z();
        new_inf.hinge_point.x() = orig_inf.hinge_point.x();
        new_inf.hinge_point.y() = -orig_inf.hinge_point.y() + y_offset;
        new_inf.hinge_point.z() = orig_inf.hinge_point.z();
        new_inf.reflection_sign = orig_inf.reflection_sign;
        new_influences.push_back(new_inf);
      }
      control_influences_.push_back(new_influences);
    }
  }

  info.num_vortices =
      static_cast<Index>(vortices_.size()) - info.first_vortex_index;
  info.num_strips =
      static_cast<Index>(strips_.size()) - info.first_strip_index;

  surface_info_.push_back(info);
}

void VortexLattice::CalculateNormals() {
  RecalculateNormals();
}

void VortexLattice::RecalculateNormals() {
  for (Index j = 0; j < static_cast<Index>(strips_.size()); ++j) {
    Strip& strip = strips_[j];

    // Get first and last vortex in strip for strip normal calculation
    if (strip.num_vortices == 0) continue;

    const HorseshoeVortex& vor_le = vortices_[strip.first_vortex_index];
    const HorseshoeVortex& vor_te =
        vortices_[strip.first_vortex_index + strip.num_vortices - 1];

    // Calculate strip direction vectors
    Vec3 dx_le = vor_le.bound_vortex_right - vor_le.bound_vortex_left;
    Vec3 dx_te = vor_te.bound_vortex_right - vor_te.bound_vortex_left;
    Vec3 dx_avg = (1.0 - spanwise_axis_fraction_) * dx_le +
                  spanwise_axis_fraction_ * dx_te;

    // Strip spanwise unit vector
    Scalar dx_mag = dx_avg.norm();
    if (dx_mag > kEpsilon) {
      strip.spanwise_unit_vector = dx_avg / dx_mag;
    }

    // Calculate normals for each vortex element in strip
    for (Index k = 0; k < strip.num_vortices; ++k) {
      HorseshoeVortex& vor = vortices_[strip.first_vortex_index + k];

      // Bound vortex direction
      Vec3 bound_dir = vor.bound_vortex_right - vor.bound_vortex_left;
      Scalar bound_len = bound_dir.norm();
      if (bound_len > kEpsilon) {
        bound_dir /= bound_len;
      }

      // Direction along chord at control point
      // Normal vector comes from incidence angle minus camber slope
      // Use Trefftz plane normal (perpendicular to span in YZ plane) for rotation
      // This matches OptVL's approach: EC = (cos, -sin*ES_y, -sin*ES_z)
      Scalar ang_cp =
          strip.twist_angle - std::atan(vor.camber_slope_at_control_point);
      Vec3 chord_dir_cp;
      chord_dir_cp.x() = std::cos(ang_cp);
      chord_dir_cp.y() = -std::sin(ang_cp) * strip.trefftz_normal_y;
      chord_dir_cp.z() = -std::sin(ang_cp) * strip.trefftz_normal_z;

      // Normal at control point = chord_dir x bound_dir
      Vec3 normal_cp = chord_dir_cp.cross(bound_dir);
      Scalar normal_mag = normal_cp.norm();
      if (normal_mag > kEpsilon) {
        vor.normal_at_control_point = normal_cp / normal_mag;
      } else {
        // Fallback to strip normal
        vor.normal_at_control_point.x() = 0.0;
        vor.normal_at_control_point.y() = strip.trefftz_normal_y;
        vor.normal_at_control_point.z() = strip.trefftz_normal_z;
      }

      // Direction along chord at vortex
      // Use Trefftz plane normal for twist rotation (same as OptVL)
      Scalar ang_vr =
          strip.twist_angle - std::atan(vor.camber_slope_at_vortex);
      Vec3 chord_dir_vr;
      chord_dir_vr.x() = std::cos(ang_vr);
      chord_dir_vr.y() = -std::sin(ang_vr) * strip.trefftz_normal_y;
      chord_dir_vr.z() = -std::sin(ang_vr) * strip.trefftz_normal_z;

      // Normal at vortex = chord_dir x bound_dir
      Vec3 normal_vr = chord_dir_vr.cross(bound_dir);
      normal_mag = normal_vr.norm();
      if (normal_mag > kEpsilon) {
        vor.normal_at_vortex = normal_vr / normal_mag;
      } else {
        vor.normal_at_vortex.x() = 0.0;
        vor.normal_at_vortex.y() = strip.trefftz_normal_y;
        vor.normal_at_vortex.z() = strip.trefftz_normal_z;
      }
    }
  }
}

void VortexLattice::CalculateChordwiseSpacing(Index num_chordwise,
                                              SpacingType spacing_type,
                                              Scalar cl_alpha_factor,
                                              std::vector<Scalar>& x_vortex,
                                              std::vector<Scalar>& x_control) const {
  x_vortex.resize(num_chordwise);
  x_control.resize(num_chordwise);

  // Generate panel edge spacing
  std::vector<Scalar> x_pt =
      GenerateCosineSpacing(num_chordwise + 1, spacing_type);

  // For VLM, the bound vortex is at the 1/4 chord of each panel,
  // and the control point is at the 3/4 chord of each panel.
  // Modified by CL-alpha factor to shift vortex position.
  for (Index i = 0; i < num_chordwise; ++i) {
    Scalar x1 = x_pt[i];
    Scalar x2 = x_pt[i + 1];
    Scalar dx = x2 - x1;

    // Standard VLM: vortex at 1/4, control at 3/4
    // With CL-alpha factor: vortex shifted aft for higher camber
    Scalar vortex_fraction = 0.25;

    // Adjust vortex position based on local CL-alpha factor
    // This accounts for higher lift-curve slope due to camber
    if (cl_alpha_factor != 1.0 && i == 0) {
      // Only adjust leading edge panel significantly
      vortex_fraction = 0.25 / cl_alpha_factor;
      vortex_fraction = std::max(0.1, std::min(0.4, vortex_fraction));
    }

    x_vortex[i] = x1 + vortex_fraction * dx;
    x_control[i] = x1 + 0.75 * dx;
  }
}

}  // namespace vlm
