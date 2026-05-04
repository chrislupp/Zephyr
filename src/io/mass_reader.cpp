/*
    Project Pegasus
    Version 0.2.0

    Copyright 2024-2026, Christopher A. Lupp

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

#include "pegasus/io/mass_reader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace pegasus {

namespace {

// Trim whitespace from both ends of a string
std::string Trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

// Remove comment (everything after ! or #)
std::string RemoveComment(const std::string& line) {
    size_t pos = line.find_first_of("!#");
    if (pos != std::string::npos) {
        return line.substr(0, pos);
    }
    return line;
}

// Extract comment text (after ! or #)
std::string ExtractComment(const std::string& line) {
    size_t pos = line.find_first_of("!#");
    if (pos != std::string::npos && pos + 1 < line.size()) {
        return Trim(line.substr(pos + 1));
    }
    return "";
}

// Convert string to lowercase
std::string ToLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

// Parse a number from string, handling various formats
bool ParseNumber(const std::string& str, Scalar& value) {
    try {
        size_t pos;
        value = std::stod(str, &pos);
        return pos == str.size();
    } catch (...) {
        return false;
    }
}

}  // namespace

// ============================================================================
// MassFileData Implementation
// ============================================================================

MassProperties MassFileData::ComputeTotalProperties() const {
    MassProperties props;

    if (components.empty()) {
        return props;
    }

    // First pass: compute total mass and CG location
    Scalar total_mass = 0.0;
    Vec3 cg_weighted_sum = Vec3::Zero();

    for (const auto& comp : components) {
        Scalar m = comp.mass * mass_unit;
        Vec3 pos = comp.position * length_unit;
        total_mass += m;
        cg_weighted_sum += m * pos;
    }

    if (total_mass <= 0.0) {
        return props;
    }

    props.mass = total_mass;
    props.cg = cg_weighted_sum / total_mass;

    // Second pass: compute inertia tensor about CG using parallel axis theorem
    // I_cg = sum(I_i + m_i * [(r_i - cg)^T(r_i - cg)*I_3x3 - (r_i - cg)(r_i -
    // cg)^T])
    Mat3 inertia = Mat3::Zero();

    Scalar length_sq = length_unit * length_unit;
    Scalar inertia_unit = mass_unit * length_sq;

    for (const auto& comp : components) {
        Scalar m = comp.mass * mass_unit;
        Vec3 pos = comp.position * length_unit;
        Vec3 r = pos - props.cg;

        // Component's own inertia tensor (in SI units)
        Mat3 I_comp;
        I_comp << comp.Ixx, -comp.Ixy, -comp.Ixz, -comp.Ixy, comp.Iyy,
            -comp.Iyz, -comp.Ixz, -comp.Iyz, comp.Izz;
        I_comp *= inertia_unit;

        // Parallel axis theorem contribution
        Scalar r_dot_r = r.squaredNorm();
        Mat3 r_outer = r * r.transpose();
        Mat3 I_parallel = m * (r_dot_r * Mat3::Identity() - r_outer);

        inertia += I_comp + I_parallel;
    }

    props.Ixx = inertia(0, 0);
    props.Iyy = inertia(1, 1);
    props.Izz = inertia(2, 2);
    props.Ixy =
        -inertia(0, 1);  // Convention: products stored with opposite sign
    props.Ixz = -inertia(0, 2);
    props.Iyz = -inertia(1, 2);

    return props;
}

Scalar MassFileData::GetTotalMass() const {
    Scalar total = 0.0;
    for (const auto& comp : components) {
        total += comp.mass * mass_unit;
    }
    return total;
}

Vec3 MassFileData::GetCenterOfGravity() const {
    Scalar total_mass = 0.0;
    Vec3 cg_weighted_sum = Vec3::Zero();

    for (const auto& comp : components) {
        Scalar m = comp.mass * mass_unit;
        Vec3 pos = comp.position * length_unit;
        total_mass += m;
        cg_weighted_sum += m * pos;
    }

    if (total_mass <= 0.0) {
        return Vec3::Zero();
    }

    return cg_weighted_sum / total_mass;
}

// ============================================================================
// MassReader Implementation
// ============================================================================

MassFileData MassReader::ReadFile(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open mass file: " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return ReadString(buffer.str());
}

MassFileData MassReader::ReadString(const std::string& content) const {
    MassFileData data;
    Vec3 scale_factors = Vec3::Ones();
    Vec3 offset_factors = Vec3::Zero();

    std::istringstream stream(content);
    std::string line;

    while (std::getline(stream, line)) {
        ParseLine(line, data, scale_factors, offset_factors);
    }

    return data;
}

void MassReader::ParseLine(const std::string& line, MassFileData& data,
                           Vec3& scale_factors, Vec3& offset_factors) const {
    std::string trimmed = Trim(line);

    // Skip empty lines and pure comment lines
    if (trimmed.empty() || trimmed[0] == '!' || trimmed[0] == '#') {
        return;
    }

    // Try parsing different line types in order
    if (ParseUnitLine(trimmed, data))
        return;
    if (ParseConstantLine(trimmed, data))
        return;
    if (ParseScaleLine(trimmed, scale_factors))
        return;
    if (ParseOffsetLine(trimmed, offset_factors))
        return;
    ParseComponentLine(line, data, scale_factors, offset_factors);
}

bool MassReader::ParseUnitLine(const std::string& line,
                               MassFileData& data) const {
    // Format: Lunit = value [unit_name]
    //         Munit = value [unit_name]
    //         Tunit = value [unit_name]

    size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos)
        return false;

    std::string key = Trim(line.substr(0, eq_pos));
    std::string value_part = Trim(RemoveComment(line.substr(eq_pos + 1)));

    // Extract numeric value (first token)
    std::istringstream iss(value_part);
    Scalar value;
    if (!(iss >> value))
        return false;

    std::string key_lower = ToLower(key);
    if (key_lower == "lunit") {
        data.length_unit = value;
        return true;
    } else if (key_lower == "munit") {
        data.mass_unit = value;
        return true;
    } else if (key_lower == "tunit") {
        data.time_unit = value;
        return true;
    }

    return false;
}

bool MassReader::ParseConstantLine(const std::string& line,
                                   MassFileData& data) const {
    // Format: g = value
    //         rho = value

    size_t eq_pos = line.find('=');
    if (eq_pos == std::string::npos)
        return false;

    std::string key = Trim(line.substr(0, eq_pos));
    std::string value_str = Trim(RemoveComment(line.substr(eq_pos + 1)));

    Scalar value;
    if (!ParseNumber(value_str, value))
        return false;

    std::string key_lower = ToLower(key);
    if (key_lower == "g") {
        data.gravity = value;
        return true;
    } else if (key_lower == "rho") {
        data.density = value;
        return true;
    }

    return false;
}

bool MassReader::ParseScaleLine(const std::string& line,
                                Vec3& scale_factors) const {
    // Format: *  sx  sy  sz  (multipliers for subsequent component positions)
    std::string trimmed = Trim(RemoveComment(line));
    if (trimmed.empty() || trimmed[0] != '*')
        return false;

    std::istringstream iss(trimmed.substr(1));
    Scalar sx, sy, sz;
    if (iss >> sx >> sy >> sz) {
        scale_factors = Vec3(sx, sy, sz);
        return true;
    }

    return false;
}

bool MassReader::ParseOffsetLine(const std::string& line,
                                 Vec3& offset_factors) const {
    // Format: +  ox  oy  oz  (offsets added to subsequent component positions)
    std::string trimmed = Trim(RemoveComment(line));
    if (trimmed.empty() || trimmed[0] != '+')
        return false;

    std::istringstream iss(trimmed.substr(1));
    Scalar ox, oy, oz;
    if (iss >> ox >> oy >> oz) {
        offset_factors = Vec3(ox, oy, oz);
        return true;
    }

    return false;
}

bool MassReader::ParseComponentLine(const std::string& line, MassFileData& data,
                                    const Vec3& scale_factors,
                                    const Vec3& offset_factors) const {
    // Format: mass  x  y  z  Ixx  Iyy  Izz  [Ixy  Ixz  Iyz]  ! comment
    // Minimum 7 values (mass, x, y, z, Ixx, Iyy, Izz), optionally 10

    std::string data_part = RemoveComment(line);
    std::string trimmed = Trim(data_part);

    if (trimmed.empty())
        return false;

    // Skip lines starting with special characters or letters (keywords)
    char first = trimmed[0];
    if (first == '*' || first == '+' || first == '!' || first == '#') {
        return false;
    }
    if (std::isalpha(static_cast<unsigned char>(first))) {
        return false;
    }

    std::istringstream iss(trimmed);
    std::vector<Scalar> values;
    Scalar v;
    while (iss >> v) {
        values.push_back(v);
    }

    // Need at least 7 values: mass, x, y, z, Ixx, Iyy, Izz
    if (values.size() < 7)
        return false;

    MassComponent comp;
    comp.mass = values[0];

    // Apply scale and offset to position
    comp.position.x() = values[1] * scale_factors.x() + offset_factors.x();
    comp.position.y() = values[2] * scale_factors.y() + offset_factors.y();
    comp.position.z() = values[3] * scale_factors.z() + offset_factors.z();

    comp.Ixx = values[4];
    comp.Iyy = values[5];
    comp.Izz = values[6];

    // Optional products of inertia
    if (values.size() >= 8)
        comp.Ixy = values[7];
    if (values.size() >= 9)
        comp.Ixz = values[8];
    if (values.size() >= 10)
        comp.Iyz = values[9];

    // Extract component name from comment
    comp.name = ExtractComment(line);

    data.components.push_back(comp);
    return true;
}

}  // namespace pegasus
