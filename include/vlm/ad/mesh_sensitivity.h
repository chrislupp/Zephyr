/*
    Project VLM
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

#pragma once

/// @file mesh_sensitivity.h
/// @brief Infrastructure for computing mesh geometry sensitivities.
///
/// This module provides the chain rule from section parameters to mesh
/// positions, matching OptVL's approach. The key relationships are:
///
/// Section Parameters → Mesh Positions:
///   RLE(strip) = (1-FC)*LE_left + FC*LE_right
///   CHORD(strip) = (1-FC)*CHORD_left + FC*CHORD_right
///   RV1(vor) = RLE1(strip) + XVR*CHORD1(strip)
///   RV2(vor) = RLE2(strip) + XVR*CHORD2(strip)
///   RC(vor) = RLE(strip) + XCP*CHORD(strip)
///
/// Where FC is the interpolation fraction within the section interval.

#include "vlm/core/types.h"
#include "vlm/utils/math_utilities.h"

#include <cmath>
#include <vector>

namespace vlm {

// Forward declarations
class VortexLattice;
class Aircraft;

namespace ad {

/// @brief Information about which strips belong to each section interval.
///
/// This matches OptVL's ICNTSEC/IPTLOC arrays. For a section interval
/// [section_idx, section_idx+1], the affected strips are:
///   [first_strip_in_interval[section_idx],
///   first_strip_in_interval[section_idx+1])
struct SectionStripMapping {
    Index surface_index{0};
    Index num_sections{0};
    Index first_strip{0};                 ///< First strip for this surface
    Index num_strips{0};                  ///< Total strips for this surface
    std::vector<Index> strip_at_section;  ///< Strip edge index at each section

    /// Get the range of strips for a section interval [isec, isec+1]
    /// @return pair of (first_strip_in_interval, num_strips_in_interval)
    [[nodiscard]] std::pair<Index, Index> GetStripRange(
        Index section_idx) const {
        if (section_idx + 1 >= static_cast<Index>(strip_at_section.size())) {
            return {0, 0};
        }
        Index first = first_strip + strip_at_section[section_idx];
        Index last = first_strip + strip_at_section[section_idx + 1];
        return {first, last - first};
    }
};

/// @brief Mesh arrays extracted for gradient computation.
///
/// These are the positions and normals that enter the AIC computation.
/// Matches OptVL's RV1, RV2, RC, ENC arrays.
struct MeshArrays {
    std::vector<double> rv1;         ///< Left bound vortex points [3*N]
    std::vector<double> rv2;         ///< Right bound vortex points [3*N]
    std::vector<double> rc;          ///< Control points [3*N]
    std::vector<double> normal;      ///< Normal vectors at control points [3*N]
    std::vector<double> chord;       ///< Local chord at each vortex [N]
    std::vector<int> surface_idx;    ///< Surface index for each vortex [N]
    std::vector<int> component_idx;  ///< Component index for each vortex [N]
    Index num_vortices;
};

/// @brief Mesh position gradients (adjoint seeds from AIC/force computation).
///
/// These accumulate d(objective)/d(mesh position).
struct MeshGradients {
    std::vector<double> d_rv1;     ///< [3*N]
    std::vector<double> d_rv2;     ///< [3*N]
    std::vector<double> d_rc;      ///< [3*N]
    std::vector<double> d_normal;  ///< [3*N]

    void Resize(Index num_vortices) {
        d_rv1.assign(3 * num_vortices, 0.0);
        d_rv2.assign(3 * num_vortices, 0.0);
        d_rc.assign(3 * num_vortices, 0.0);
        d_normal.assign(3 * num_vortices, 0.0);
    }

    void Zero() {
        std::fill(d_rv1.begin(), d_rv1.end(), 0.0);
        std::fill(d_rv2.begin(), d_rv2.end(), 0.0);
        std::fill(d_rc.begin(), d_rc.end(), 0.0);
        std::fill(d_normal.begin(), d_normal.end(), 0.0);
    }
};

/// @brief Section parameter gradients (output of chain rule).
///
/// These are the final gradients w.r.t. geometry parameters.
/// Matches OptVL's XYZLES_diff, CHORDS_diff, AINCS_diff.
struct SectionGradients {
    double d_le_x{0};   ///< d/d(leading edge x)
    double d_le_y{0};   ///< d/d(leading edge y)
    double d_le_z{0};   ///< d/d(leading edge z)
    double d_chord{0};  ///< d/d(chord)
    double d_twist{0};  ///< d/d(twist angle)

    SectionGradients() = default;
};

// ============================================================================
// Interpolation Formulas (Forward Mode)
// ============================================================================

/// @brief Compute strip interpolation fraction for a given strip.
///
/// Given a strip within a section interval, compute the fraction FC
/// where 0 = left section, 1 = right section.
/// Uses the mesh generator's spacing to get correct fractions.
///
/// @param strip_idx_in_interval Strip index within the section interval
/// (0-based)
/// @param num_strips_in_interval Total strips in this section interval
/// @param spacing Spacing type used by the mesh generator
/// @return Interpolation fraction FC in [0, 1]
inline double ComputeStripInterpolationFraction(
    Index strip_idx_in_interval, Index num_strips_in_interval,
    SpacingType spacing = SpacingType::kCosine) {
    if (num_strips_in_interval <= 1)
        return 0.5;
    // The mesh generator uses GenerateCosineSpacing(2*N+1, spacing)
    // where odd indices [2n+1] are control point (center) fractions
    std::vector<Scalar> fracs =
        GenerateCosineSpacing(2 * num_strips_in_interval + 1, spacing);
    return fracs[2 * strip_idx_in_interval + 1];
}

/// @brief Compute edge interpolation fractions for a strip.
///
/// Returns F1 (left edge) and F2 (right edge) fractions.
/// Uses the mesh generator's spacing to get correct fractions.
inline void ComputeStripEdgeFractions(
    Index strip_idx_in_interval, Index num_strips_in_interval, double* f1,
    double* f2, SpacingType spacing = SpacingType::kCosine) {
    // The mesh generator uses GenerateCosineSpacing(2*N+1, spacing)
    // where even indices [2n] are panel edge fractions
    std::vector<Scalar> fracs =
        GenerateCosineSpacing(2 * num_strips_in_interval + 1, spacing);
    *f1 = fracs[2 * strip_idx_in_interval];
    *f2 = fracs[2 * strip_idx_in_interval + 2];
}

/// @brief Pre-compute all spanwise fractions for an interval.
///
/// More efficient than calling ComputeStripEdgeFractions per strip.
/// Returns arrays of f1[N], f2[N], fc[N] for N strips.
inline void BuildSpanwiseFractions(Index num_strips, SpacingType spacing,
                                   std::vector<double>& f1_out,
                                   std::vector<double>& f2_out,
                                   std::vector<double>& fc_out) {
    std::vector<Scalar> fracs =
        GenerateCosineSpacing(2 * num_strips + 1, spacing);
    f1_out.resize(num_strips);
    f2_out.resize(num_strips);
    fc_out.resize(num_strips);
    for (Index i = 0; i < num_strips; ++i) {
        f1_out[i] = fracs[2 * i];
        f2_out[i] = fracs[2 * i + 1];  // Note: fc is at odd index (center)
        fc_out[i] = fracs[2 * i + 1];
    }
    // f2 uses the next even index
    for (Index i = 0; i < num_strips; ++i) {
        f2_out[i] = fracs[2 * i + 2];
    }
}

// ============================================================================
// Gradient Chain Rule (Reverse Mode)
// ============================================================================

/// @brief Propagate mesh gradients back to section leading edge.
///
/// Given d_rle (gradient of objective w.r.t. strip leading edge),
/// accumulate contributions to left and right section LE gradients.
///
/// Chain rule: RLE = (1-FC)*LE_left + FC*LE_right
/// Adjoint:    d_LE_left += (1-FC)*d_RLE
///             d_LE_right += FC*d_RLE
///
/// @param d_rle Gradient w.r.t. strip leading edge [3]
/// @param fc Interpolation fraction
/// @param d_le_left Output: accumulated gradient for left section LE [3]
/// @param d_le_right Output: accumulated gradient for right section LE [3]
inline void PropagateLeadingEdgeGradient(const double* d_rle, double fc,
                                         double* d_le_left,
                                         double* d_le_right) {
    for (int i = 0; i < 3; ++i) {
        d_le_left[i] += (1.0 - fc) * d_rle[i];
        d_le_right[i] += fc * d_rle[i];
    }
}

/// @brief Propagate mesh gradients back to section chord.
///
/// Chain rule: CHORD(strip) = (1-FC)*CHORD_left + FC*CHORD_right
///             RV1.x = RLE1.x + XVR*CHORD1
///             RC.x = RLE.x + XCP*CHORD
///
/// Adjoint: d_CHORD_left += (1-F1)*XVR*d_RV1.x + (1-FC)*XCP*d_RC.x
///          d_CHORD_right += F1*XVR*d_RV1.x + FC*XCP*d_RC.x
///
/// @param d_rv1 Gradient w.r.t. left bound vortex [3]
/// @param d_rv2 Gradient w.r.t. right bound vortex [3]
/// @param d_rc Gradient w.r.t. control point [3]
/// @param xvr Chordwise vortex position fraction
/// @param xcp Chordwise control point fraction
/// @param f1 Left edge interpolation fraction
/// @param f2 Right edge interpolation fraction
/// @param fc Center interpolation fraction
/// @param d_chord_left Output: accumulated chord gradient for left section
/// @param d_chord_right Output: accumulated chord gradient for right section
inline void PropagateChordGradient(const double* d_rv1, const double* d_rv2,
                                   const double* d_rc, double xvr, double xcp,
                                   double f1, double f2, double fc,
                                   double* d_chord_left,
                                   double* d_chord_right) {
    // Contribution from RV1 (left edge)
    *d_chord_left += (1.0 - f1) * xvr * d_rv1[0];
    *d_chord_right += f1 * xvr * d_rv1[0];

    // Contribution from RV2 (right edge)
    *d_chord_left += (1.0 - f2) * xvr * d_rv2[0];
    *d_chord_right += f2 * xvr * d_rv2[0];

    // Contribution from RC (center)
    *d_chord_left += (1.0 - fc) * xcp * d_rc[0];
    *d_chord_right += fc * xcp * d_rc[0];
}

/// @brief Propagate normal vector gradients back to section twist.
///
/// The normal vector depends on twist through:
///   angle = twist - atan(camber_slope)
///   EC = (cos(angle), -sin(angle)*ES_y, -sin(angle)*ES_z)
///   ENC = normalize(EC × EB)
///
/// This is a complex function, so we compute it analytically here.
///
/// @param d_normal Gradient w.r.t. normal vector [3]
/// @param twist Current twist angle (radians)
/// @param trefftz_normal_y Y component of Trefftz plane normal
/// @param trefftz_normal_z Z component of Trefftz plane normal
/// @param bound_dir Bound vortex direction [3]
/// @param camber_slope Camber slope at this location
/// @return Gradient w.r.t. twist angle
inline double PropagateTwistGradient(const double* d_normal, double twist,
                                     double trefftz_normal_y,
                                     double trefftz_normal_z,
                                     const double* bound_dir,
                                     double camber_slope) {
    // Effective angle
    double angle = twist - std::atan(camber_slope);
    double ca = std::cos(angle);
    double sa = std::sin(angle);

    // Chord direction EC
    double ec[3];
    ec[0] = ca;
    ec[1] = -sa * trefftz_normal_y;
    ec[2] = -sa * trefftz_normal_z;

    // d(EC)/d(angle)
    double dec_dangle[3];
    dec_dangle[0] = -sa;
    dec_dangle[1] = -ca * trefftz_normal_y;
    dec_dangle[2] = -ca * trefftz_normal_z;

    // Normal = normalize(EC × EB)
    // Let ECXB = EC × EB
    double ecxb[3];
    ecxb[0] = ec[1] * bound_dir[2] - ec[2] * bound_dir[1];
    ecxb[1] = ec[2] * bound_dir[0] - ec[0] * bound_dir[2];
    ecxb[2] = ec[0] * bound_dir[1] - ec[1] * bound_dir[0];

    double emag =
        std::sqrt(ecxb[0] * ecxb[0] + ecxb[1] * ecxb[1] + ecxb[2] * ecxb[2]);
    if (emag < 1e-12)
        return 0.0;

    // Normal = ECXB / emag
    double normal[3];
    normal[0] = ecxb[0] / emag;
    normal[1] = ecxb[1] / emag;
    normal[2] = ecxb[2] / emag;

    // d(ECXB)/d(angle) = d(EC)/d(angle) × EB
    double decxb_dangle[3];
    decxb_dangle[0] =
        dec_dangle[1] * bound_dir[2] - dec_dangle[2] * bound_dir[1];
    decxb_dangle[1] =
        dec_dangle[2] * bound_dir[0] - dec_dangle[0] * bound_dir[2];
    decxb_dangle[2] =
        dec_dangle[0] * bound_dir[1] - dec_dangle[1] * bound_dir[0];

    // d(emag)/d(angle) = normal · d(ECXB)/d(angle)
    double demag_dangle = normal[0] * decxb_dangle[0] +
                          normal[1] * decxb_dangle[1] +
                          normal[2] * decxb_dangle[2];

    // d(normal)/d(angle) = (d(ECXB)/d(angle) - normal * d(emag)/d(angle)) /
    // emag
    double dnormal_dangle[3];
    dnormal_dangle[0] = (decxb_dangle[0] - normal[0] * demag_dangle) / emag;
    dnormal_dangle[1] = (decxb_dangle[1] - normal[1] * demag_dangle) / emag;
    dnormal_dangle[2] = (decxb_dangle[2] - normal[2] * demag_dangle) / emag;

    // Chain rule: d(obj)/d(twist) = d(obj)/d(normal) · d(normal)/d(twist)
    // Note: d(angle)/d(twist) = 1
    return d_normal[0] * dnormal_dangle[0] + d_normal[1] * dnormal_dangle[1] +
           d_normal[2] * dnormal_dangle[2];
}

/// @brief Complete section gradient computation for one section.
///
/// This function takes mesh gradients (from AIC/force adjoint) and
/// propagates them back to section parameters using the chain rule.
///
/// @param mesh_grads Input: d(objective)/d(mesh positions)
/// @param mapping Section-to-strip mapping for this surface
/// @param section_idx Which section to compute gradients for
/// @param strip_data Additional strip data (chords, twists, etc.)
/// @param chordwise_fractions Chordwise position fractions (XVR, XCP)
/// @return Gradients w.r.t. this section's parameters
SectionGradients ComputeSectionGradientFromMesh(
    const MeshGradients& mesh_grads, const SectionStripMapping& mapping,
    Index section_idx, const VortexLattice& mesh,
    const std::vector<double>& xvr_fractions,
    const std::vector<double>& xcp_fractions);

// ============================================================================
// Mapping Construction
// ============================================================================

/// @brief Build section-to-strip mapping for a surface.
///
/// This reconstructs the mapping from Aircraft + VortexLattice data.
/// Matches OptVL's ICNTSEC/IPTLOC logic.
///
/// @param aircraft The aircraft geometry
/// @param mesh The vortex lattice mesh
/// @param surface_idx Surface index
/// @return The section-to-strip mapping
SectionStripMapping BuildSectionStripMapping(const Aircraft& aircraft,
                                             const VortexLattice& mesh,
                                             Index surface_idx);

/// @brief Extract mesh arrays for gradient computation.
///
/// Extracts raw position and normal arrays from VortexLattice into
/// contiguous arrays suitable for analytical differentiation.
///
/// @param mesh The vortex lattice mesh
/// @return Extracted mesh arrays
MeshArrays ExtractMeshArrays(const VortexLattice& mesh);

/// @brief Build default chordwise fractions.
///
/// Returns the standard VLM fractions: vortex at 1/4 panel, control at 3/4.
///
/// @param num_chordwise Number of chordwise panels
/// @param xvr_out Output vortex position fractions
/// @param xcp_out Output control point fractions
void BuildDefaultChordwiseFractions(Index num_chordwise,
                                    std::vector<double>& xvr_out,
                                    std::vector<double>& xcp_out,
                                    SpacingType spacing = SpacingType::kCosine);

}  // namespace ad
}  // namespace vlm
