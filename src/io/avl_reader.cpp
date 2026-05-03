// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/io/avl_reader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace vlm {

namespace {

/// @brief Trims whitespace from both ends of a string.
std::string Trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
        return "";
    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

/// @brief Checks if a line is a comment or empty.
bool IsCommentOrEmpty(const std::string& line) {
    std::string trimmed = Trim(line);
    return trimmed.empty() || trimmed[0] == '#' || trimmed[0] == '!';
}

/// @brief Extracts keyword from a line (first 4 characters, uppercase).
std::string GetKeyword(const std::string& line) {
    std::string trimmed = Trim(line);
    if (trimmed.size() < 4)
        return trimmed;
    std::string keyword = trimmed.substr(0, 4);
    std::transform(keyword.begin(), keyword.end(), keyword.begin(), ::toupper);
    return keyword;
}

/// @brief Parses space/comma separated numbers from a line.
std::vector<Scalar> ParseNumbers(const std::string& line) {
    std::vector<Scalar> numbers;
    std::string cleaned = line;
    // Replace commas with spaces
    std::replace(cleaned.begin(), cleaned.end(), ',', ' ');

    std::istringstream iss(cleaned);
    Scalar value;
    while (iss >> value) {
        numbers.push_back(value);
    }
    return numbers;
}

/// @brief Gets the directory part of a file path.
std::string GetDirectory(const std::string& filepath) {
    size_t pos = filepath.find_last_of("/\\");
    if (pos == std::string::npos)
        return ".";
    return filepath.substr(0, pos);
}

/// @brief Resolves a relative path against a base directory.
std::string ResolvePath(const std::string& base_path,
                        const std::string& relative_path) {
    if (relative_path.empty())
        return relative_path;
    // If absolute path, return as-is
    if (relative_path[0] == '/' ||
        (relative_path.size() > 1 && relative_path[1] == ':')) {
        return relative_path;
    }
    if (base_path.empty())
        return relative_path;
    return base_path + "/" + relative_path;
}

/// @brief Reads lines from a stream, skipping comments.
class LineReader {
public:
    explicit LineReader(std::istream& stream) : stream_(stream) {}

    bool ReadLine(std::string& line) {
        while (std::getline(stream_, line)) {
            ++line_number_;
            if (!IsCommentOrEmpty(line)) {
                return true;
            }
        }
        return false;
    }

    int GetLineNumber() const { return line_number_; }

private:
    std::istream& stream_;
    int line_number_{0};
};

}  // namespace

Aircraft AvlReader::ReadFile(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("AvlReader::ReadFile: cannot open file '" +
                                 filename + "'");
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return ParseString(buffer.str(), GetDirectory(filename));
}

Aircraft AvlReader::ParseString(const std::string& content,
                                const std::string& base_path) const {
    std::istringstream stream(content);
    LineReader reader(stream);
    std::string line;

    Aircraft aircraft;

    // Line 1: Title
    if (!reader.ReadLine(line)) {
        throw std::runtime_error("AvlReader: unexpected end of file (title)");
    }
    aircraft.SetTitle(Trim(line));
    if (verbose_) {
        std::cout << "Configuration: " << aircraft.GetTitle() << std::endl;
    }

    // Line 2: Mach number
    if (!reader.ReadLine(line)) {
        throw std::runtime_error("AvlReader: unexpected end of file (Mach)");
    }
    auto mach_values = ParseNumbers(line);
    if (!mach_values.empty()) {
        aircraft.SetMachNumber(mach_values[0]);
    }

    // Line 3: Symmetry (IYsym, IZsym, Zsym)
    if (!reader.ReadLine(line)) {
        throw std::runtime_error(
            "AvlReader: unexpected end of file (symmetry)");
    }
    auto sym_values = ParseNumbers(line);
    if (sym_values.size() >= 3) {
        int iy_sym = static_cast<int>(sym_values[0]);
        int iz_sym = static_cast<int>(sym_values[1]);
        Scalar z_sym = sym_values[2];

        if (iy_sym > 0)
            aircraft.SetYSymmetry(SymmetryType::kSymmetric);
        else if (iy_sym < 0)
            aircraft.SetYSymmetry(SymmetryType::kAntisymmetric);

        if (iz_sym > 0)
            aircraft.SetZSymmetry(SymmetryType::kSymmetric);
        else if (iz_sym < 0)
            aircraft.SetZSymmetry(SymmetryType::kAntisymmetric);

        aircraft.SetZSymmetryPlane(z_sym);
    }

    // Line 4: Reference quantities (Sref, Cref, Bref)
    if (!reader.ReadLine(line)) {
        throw std::runtime_error(
            "AvlReader: unexpected end of file (reference)");
    }
    auto ref_values = ParseNumbers(line);
    if (ref_values.size() >= 3) {
        aircraft.SetReferenceArea(ref_values[0]);
        aircraft.SetReferenceChord(ref_values[1]);
        aircraft.SetReferenceSpan(ref_values[2]);
    }

    // Line 5: Moment reference point (Xref, Yref, Zref)
    if (!reader.ReadLine(line)) {
        throw std::runtime_error(
            "AvlReader: unexpected end of file (moment ref)");
    }
    auto mom_ref_values = ParseNumbers(line);
    if (mom_ref_values.size() >= 3) {
        aircraft.SetMomentReferencePoint(
            Vec3(mom_ref_values[0], mom_ref_values[1], mom_ref_values[2]));
    }

    // Line 6: Profile drag coefficient (optional)
    if (reader.ReadLine(line)) {
        std::string keyword = GetKeyword(line);
        if (keyword != "SURF" && keyword != "BODY") {
            auto cdp_values = ParseNumbers(line);
            if (!cdp_values.empty()) {
                aircraft.SetProfileDragCoefficient(cdp_values[0]);
            }
        } else {
            // This line is actually a SURFACE or BODY keyword, process it
            goto process_keyword;
        }
    }

    // Parse surfaces and bodies
    while (reader.ReadLine(line)) {
    process_keyword:
        std::string keyword = GetKeyword(line);

        if (keyword == "SURF") {
            // SURFACE block
            Surface surface;

            // Surface name is on the next line
            if (!reader.ReadLine(line)) {
                throw std::runtime_error(
                    "AvlReader: unexpected end of file after SURFACE");
            }
            surface.SetName(Trim(line));
            if (verbose_) {
                std::cout << "  Surface: " << surface.GetName() << std::endl;
            }

            // Paneling parameters (Nchord, Cspace, [Nspan, Sspace])
            if (!reader.ReadLine(line)) {
                throw std::runtime_error(
                    "AvlReader: unexpected end of file (paneling)");
            }
            auto panel_values = ParseNumbers(line);
            if (panel_values.size() >= 2) {
                surface.SetChordwisePanelCount(
                    static_cast<int>(panel_values[0]));
                surface.SetChordwiseSpacing(panel_values[1]);
                if (panel_values.size() >= 4) {
                    surface.SetSpanwisePanelCount(
                        static_cast<int>(panel_values[2]));
                    surface.SetSpanwiseSpacing(panel_values[3]);
                }
            }

            // Parse surface modifiers and sections
            Section current_section;
            bool has_section = false;

            while (reader.ReadLine(line)) {
                keyword = GetKeyword(line);

                if (keyword == "SURF" || keyword == "BODY") {
                    // New surface or body - save current surface and reprocess
                    if (has_section) {
                        surface.AddSection(current_section);
                    }
                    aircraft.AddSurface(std::move(surface));
                    goto process_keyword;
                } else if (keyword == "SCAL") {
                    // SCALE
                    if (reader.ReadLine(line)) {
                        auto scale = ParseNumbers(line);
                        if (scale.size() >= 3) {
                            surface.SetScale(
                                Vec3(scale[0], scale[1], scale[2]));
                        }
                    }
                } else if (keyword == "TRAN") {
                    // TRANSLATE
                    if (reader.ReadLine(line)) {
                        auto trans = ParseNumbers(line);
                        if (trans.size() >= 3) {
                            surface.SetTranslation(
                                Vec3(trans[0], trans[1], trans[2]));
                        }
                    }
                } else if (keyword == "ANGL") {
                    // ANGLE
                    if (reader.ReadLine(line)) {
                        auto angle = ParseNumbers(line);
                        if (!angle.empty()) {
                            surface.SetAngleOffset(angle[0]);
                        }
                    }
                } else if (keyword == "YDUP") {
                    // YDUPLICATE
                    if (reader.ReadLine(line)) {
                        auto ydup = ParseNumbers(line);
                        if (!ydup.empty()) {
                            surface.SetYDuplicate(ydup[0]);
                        }
                    }
                } else if (keyword == "COMP") {
                    // COMPONENT
                    if (reader.ReadLine(line)) {
                        auto comp = ParseNumbers(line);
                        if (!comp.empty()) {
                            surface.SetComponentIndex(
                                static_cast<int>(comp[0]));
                        }
                    }
                } else if (keyword == "NOWA") {
                    // NOWAKE
                    surface.SetHasWake(false);
                } else if (keyword == "NOLO") {
                    // NOLOAD
                    surface.SetContributesToLoads(false);
                } else if (keyword == "NOAL") {
                    // NOALBE
                    surface.SetAffectedByAlphaBeta(false);
                } else if (keyword == "CDCL") {
                    // CDCL - drag polar
                    if (reader.ReadLine(line)) {
                        auto polar = ParseNumbers(line);
                        std::array<Scalar, 6> polar_array = {0, 0, 0, 0, 0, 0};
                        for (size_t i = 0;
                             i < std::min(polar.size(), size_t(6)); ++i) {
                            polar_array[i] = polar[i];
                        }
                        surface.SetClCdPolar(polar_array);
                    }
                } else if (keyword == "SECT") {
                    // SECTION
                    if (has_section) {
                        surface.AddSection(current_section);
                    }
                    current_section = Section();
                    has_section = true;

                    // Section data (Xle, Yle, Zle, Chord, Ainc, [Nspan,
                    // Sspace])
                    if (reader.ReadLine(line)) {
                        auto sec_data = ParseNumbers(line);
                        if (sec_data.size() >= 5) {
                            current_section.SetLeadingEdgePosition(
                                Vec3(sec_data[0], sec_data[1], sec_data[2]));
                            current_section.SetChord(sec_data[3]);
                            current_section.SetTwistAngle(sec_data[4]);
                            if (sec_data.size() >= 7) {
                                current_section.SetSpanwisePanelCount(
                                    static_cast<int>(sec_data[5]));
                                current_section.SetSpanwiseSpacing(sec_data[6]);
                            }
                        }
                    }
                } else if (keyword == "NACA") {
                    // NACA airfoil
                    if (reader.ReadLine(line)) {
                        int naca_code = static_cast<int>(ParseNumbers(line)[0]);
                        auto airfoil = std::make_shared<Airfoil>(
                            Airfoil::CreateFromNaca4Digit(naca_code));
                        current_section.SetAirfoil(airfoil);
                    }
                } else if (keyword == "AFIL" || keyword == "AIRF") {
                    // AFILE - airfoil from file
                    if (reader.ReadLine(line)) {
                        std::string airfoil_path =
                            ResolvePath(base_path, Trim(line));
                        try {
                            auto airfoil = std::make_shared<Airfoil>(
                                Airfoil::CreateFromCoordinateFile(
                                    airfoil_path));
                            current_section.SetAirfoil(airfoil);
                        } catch (const std::exception& e) {
                            if (verbose_) {
                                std::cerr
                                    << "Warning: could not load airfoil file: "
                                    << airfoil_path << std::endl;
                            }
                        }
                    }
                } else if (keyword == "CLAF") {
                    // CLAF - CL-alpha factor
                    if (reader.ReadLine(line)) {
                        auto claf = ParseNumbers(line);
                        if (!claf.empty()) {
                            current_section.SetClAlphaFactor(claf[0]);
                        }
                    }
                } else if (keyword == "CONT") {
                    // CONTROL
                    if (reader.ReadLine(line)) {
                        // Parse: name gain xhinge hx hy hz SgnDup
                        std::istringstream iss(line);
                        ControlSurfaceDefinition control;
                        iss >> control.name >> control.gain >>
                            control.hinge_x_fraction;
                        Scalar hx, hy, hz;
                        iss >> hx >> hy >> hz;
                        control.hinge_vector = Vec3(hx, hy, hz);
                        iss >> control.sign;
                        current_section.AddControlSurface(control);
                    }
                } else if (keyword == "DESI") {
                    // DESIGN variable
                    if (reader.ReadLine(line)) {
                        std::istringstream iss(line);
                        DesignVariableDefinition design;
                        iss >> design.name >> design.weight;
                        current_section.AddDesignVariable(design);
                    }
                }
            }

            // End of file - save current surface
            if (has_section) {
                surface.AddSection(current_section);
            }
            aircraft.AddSurface(std::move(surface));

        } else if (keyword == "BODY") {
            // BODY block
            Body body;

            // Body name
            if (!reader.ReadLine(line)) {
                throw std::runtime_error(
                    "AvlReader: unexpected end of file after BODY");
            }
            body.SetName(Trim(line));
            if (verbose_) {
                std::cout << "  Body: " << body.GetName() << std::endl;
            }

            // Paneling parameters
            if (!reader.ReadLine(line)) {
                throw std::runtime_error(
                    "AvlReader: unexpected end of file (body paneling)");
            }
            auto panel_values = ParseNumbers(line);
            if (panel_values.size() >= 2) {
                body.SetNodeCount(static_cast<int>(panel_values[0]));
                body.SetNodeSpacing(panel_values[1]);
            }

            // Parse body modifiers and sections
            while (reader.ReadLine(line)) {
                keyword = GetKeyword(line);

                if (keyword == "SURF" || keyword == "BODY") {
                    // New surface or body
                    aircraft.AddBody(std::move(body));
                    goto process_keyword;
                } else if (keyword == "SCAL") {
                    if (reader.ReadLine(line)) {
                        auto scale = ParseNumbers(line);
                        if (scale.size() >= 3) {
                            body.SetScale(Vec3(scale[0], scale[1], scale[2]));
                        }
                    }
                } else if (keyword == "TRAN") {
                    if (reader.ReadLine(line)) {
                        auto trans = ParseNumbers(line);
                        if (trans.size() >= 3) {
                            body.SetTranslation(
                                Vec3(trans[0], trans[1], trans[2]));
                        }
                    }
                } else if (keyword == "YDUP") {
                    body.SetYDuplicate(true);
                } else if (keyword == "BSEC") {
                    // BSEC - body section
                    if (reader.ReadLine(line)) {
                        auto sec_data = ParseNumbers(line);
                        if (sec_data.size() >= 5) {
                            BodySection section;
                            section.center =
                                Vec3(sec_data[0], sec_data[1], sec_data[2]);
                            section.width = sec_data[3];
                            section.height = sec_data[4];
                            body.AddSection(section);
                        }
                    }
                }
            }

            // End of file
            aircraft.AddBody(std::move(body));
        }
    }

    return aircraft;
}

}  // namespace vlm
