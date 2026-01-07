// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#ifndef VLM_IO_AVL_READER_H_
#define VLM_IO_AVL_READER_H_

#include <string>

#include "vlm/geometry/aircraft.h"

namespace vlm {

/// @brief Reader for AVL geometry files (.avl).
///
/// Parses the Athena Vortex Lattice (AVL) geometry file format, which
/// defines aircraft configurations with surfaces, sections, and bodies.
///
/// @par File Format Overview:
/// An AVL file contains:
/// - Title line
/// - Mach number
/// - Symmetry settings (IYsym, IZsym, Zsym)
/// - Reference quantities (Sref, Cref, Bref)
/// - Moment reference point (Xref, Yref, Zref)
/// - Profile drag (CDp)
/// - SURFACE blocks containing SECTION definitions
/// - Optional BODY blocks
///
/// @par Example Usage:
/// @code
/// AvlReader reader;
/// Aircraft aircraft = reader.ReadFile("aircraft.avl");
/// @endcode
///
/// Reference: AVL User Guide by Mark Drela & Harold Youngren
class AvlReader {
 public:
  AvlReader() = default;

  /// @brief Reads an AVL file and returns an Aircraft object.
  ///
  /// @param filename Path to the AVL file
  /// @return Parsed Aircraft configuration
  /// @throws std::runtime_error if file cannot be opened or parsed
  Aircraft ReadFile(const std::string& filename) const;

  /// @brief Parses AVL content from a string.
  ///
  /// @param content AVL file content as a string
  /// @param base_path Base path for resolving relative file references
  /// @return Parsed Aircraft configuration
  Aircraft ParseString(const std::string& content,
                       const std::string& base_path = "") const;

  /// @brief Sets whether to print verbose parsing information.
  void SetVerbose(bool verbose) { verbose_ = verbose; }

 private:
  bool verbose_{false};
};

}  // namespace vlm

#endif  // VLM_IO_AVL_READER_H_
