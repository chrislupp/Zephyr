// Copyright (c) 2024 VLM Authors
// SPDX-License-Identifier: AGPL-3.0-or-later

#include "vlm/geometry/airfoil.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace vlm {

Airfoil Airfoil::CreateFlatPlate() {
  Airfoil airfoil;
  airfoil.name_ = "Flat Plate";
  airfoil.camber_x_ = {0.0, 1.0};
  airfoil.camber_y_ = {0.0, 0.0};
  airfoil.thickness_ = {0.0, 0.0};
  airfoil.max_camber_ = 0.0;
  airfoil.max_camber_position_ = 0.0;
  airfoil.max_thickness_ = 0.0;
  airfoil.BuildCamberSpline();
  return airfoil;
}

Airfoil Airfoil::CreateFromNaca4Digit(const std::string& naca_code,
                                      int num_points) {
  // Parse string to integer
  int code = 0;
  try {
    code = std::stoi(naca_code);
  } catch (const std::exception&) {
    throw std::invalid_argument(
        "Airfoil::CreateFromNaca4Digit: invalid NACA code '" + naca_code + "'");
  }
  return CreateFromNaca4Digit(code, num_points);
}

Airfoil Airfoil::CreateFromNaca4Digit(int naca_code, int num_points) {
  if (naca_code < 0 || naca_code > 9999) {
    throw std::invalid_argument(
        "Airfoil::CreateFromNaca4Digit: NACA code must be 0-9999");
  }

  // Parse NACA 4-digit code: MPTT
  // M = max camber (% chord)
  // P = position of max camber (tenths of chord)
  // TT = max thickness (% chord)
  int m_digit = naca_code / 1000;
  int p_digit = (naca_code - 1000 * m_digit) / 100;
  int tt_digits = naca_code - 1000 * m_digit - 100 * p_digit;

  Scalar max_camber = static_cast<Scalar>(m_digit) / 100.0;
  Scalar camber_position = static_cast<Scalar>(p_digit) / 10.0;
  Scalar max_thickness = static_cast<Scalar>(tt_digits) / 100.0;

  Airfoil airfoil;
  airfoil.name_ = "NACA " + std::to_string(naca_code);
  airfoil.max_camber_ = max_camber;
  airfoil.max_camber_position_ = camber_position;
  airfoil.max_thickness_ = max_thickness;

  // Generate camber line coordinates
  airfoil.camber_x_.resize(num_points);
  airfoil.camber_y_.resize(num_points);
  airfoil.thickness_.resize(num_points);

  for (int i = 0; i < num_points; ++i) {
    Scalar x = static_cast<Scalar>(i) / static_cast<Scalar>(num_points - 1);
    airfoil.camber_x_[i] = x;

    // NACA 4-digit camber line equation
    Scalar y_c = 0.0;
    if (max_camber > 0.0 && camber_position > 0.0) {
      if (x < camber_position) {
        // Forward of max camber
        y_c = max_camber / (camber_position * camber_position) *
              (2.0 * camber_position * x - x * x);
      } else {
        // Aft of max camber
        Scalar one_minus_p = 1.0 - camber_position;
        y_c = max_camber / (one_minus_p * one_minus_p) *
              ((1.0 - 2.0 * camber_position) + 2.0 * camber_position * x -
               x * x);
      }
    }
    airfoil.camber_y_[i] = y_c;

    // NACA 4-digit thickness distribution (symmetric)
    // t(x) = 5*t_max * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 +
    //                   0.2843*x^3 - 0.1015*x^4)
    // Note: original formula has -0.1036*x^4 for open trailing edge
    Scalar t = 5.0 * max_thickness *
               (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x +
                0.2843 * x * x * x - 0.1015 * x * x * x * x);
    airfoil.thickness_[i] = t;
  }

  airfoil.BuildCamberSpline();
  return airfoil;
}

Airfoil Airfoil::CreateFromCoordinateFile(const std::string& filename,
                                          int /*num_camber_points*/) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error(
        "Airfoil::CreateFromCoordinateFile: cannot open file '" + filename +
        "'");
  }

  std::vector<Scalar> x_coords;
  std::vector<Scalar> y_coords;
  std::string line;
  std::string airfoil_name;
  bool first_line = true;

  // Regex to detect lines with letters (header lines)
  std::regex alpha_regex("[a-zA-Z]");

  while (std::getline(file, line)) {
    // Skip empty lines
    if (line.empty() || std::all_of(line.begin(), line.end(), ::isspace)) {
      continue;
    }

    // Check for header line (contains letters)
    if (std::regex_search(line, alpha_regex)) {
      if (first_line) {
        airfoil_name = line;
        // Trim whitespace
        airfoil_name.erase(0, airfoil_name.find_first_not_of(" \t\r\n"));
        airfoil_name.erase(airfoil_name.find_last_not_of(" \t\r\n") + 1);
      }
      first_line = false;
      continue;
    }
    first_line = false;

    // Parse coordinate line
    std::istringstream iss(line);
    Scalar x, y;
    if (iss >> x >> y) {
      x_coords.push_back(x);
      y_coords.push_back(y);
    }
  }

  if (x_coords.size() < 3) {
    throw std::runtime_error(
        "Airfoil::CreateFromCoordinateFile: insufficient coordinates in '" +
        filename + "'");
  }

  Airfoil airfoil =
      CreateFromCoordinates(x_coords, y_coords, false);
  if (!airfoil_name.empty()) {
    airfoil.name_ = airfoil_name;
  } else {
    airfoil.name_ = filename;
  }

  return airfoil;
}

Airfoil Airfoil::CreateFromCoordinates(const std::vector<Scalar>& x_coords,
                                       const std::vector<Scalar>& y_coords,
                                       bool is_camber_line) {
  if (x_coords.size() != y_coords.size()) {
    throw std::invalid_argument(
        "Airfoil::CreateFromCoordinates: x and y must have same size");
  }
  if (x_coords.size() < 2) {
    throw std::invalid_argument(
        "Airfoil::CreateFromCoordinates: need at least 2 points");
  }

  Airfoil airfoil;

  if (is_camber_line) {
    // Direct camber line specification
    airfoil.camber_x_ = x_coords;
    airfoil.camber_y_ = y_coords;
    airfoil.thickness_.resize(x_coords.size(), 0.0);
  } else {
    // Full airfoil coordinates - extract camber line
    airfoil.ExtractCamberFromCoordinates(x_coords, y_coords);
  }

  // Compute max camber and position
  airfoil.max_camber_ = 0.0;
  airfoil.max_camber_position_ = 0.0;
  for (size_t i = 0; i < airfoil.camber_y_.size(); ++i) {
    if (std::abs(airfoil.camber_y_[i]) > std::abs(airfoil.max_camber_)) {
      airfoil.max_camber_ = airfoil.camber_y_[i];
      airfoil.max_camber_position_ = airfoil.camber_x_[i];
    }
  }

  // Compute max thickness
  airfoil.max_thickness_ = 0.0;
  for (const auto& t : airfoil.thickness_) {
    if (t > airfoil.max_thickness_) {
      airfoil.max_thickness_ = t;
    }
  }

  airfoil.BuildCamberSpline();
  return airfoil;
}

Scalar Airfoil::EvaluateCamber(Scalar x) const {
  if (!IsValid()) {
    return 0.0;
  }
  // Clamp x to [0, 1]
  x = std::max(0.0, std::min(1.0, x));
  return camber_spline_.Evaluate(x);
}

Scalar Airfoil::EvaluateCamberSlope(Scalar x) const {
  if (!IsValid()) {
    return 0.0;
  }
  // Clamp x to [0, 1]
  x = std::max(0.0, std::min(1.0, x));
  return camber_spline_.EvaluateDerivative(x);
}

Scalar Airfoil::EvaluateThickness(Scalar x) const {
  if (!thickness_spline_.IsValid()) {
    return 0.0;
  }
  x = std::max(0.0, std::min(1.0, x));
  return thickness_spline_.Evaluate(x);
}

void Airfoil::BuildCamberSpline() {
  if (camber_x_.size() >= 2) {
    camber_spline_.Build(camber_x_, camber_y_);
  }
  if (thickness_.size() >= 2 && camber_x_.size() == thickness_.size()) {
    thickness_spline_.Build(camber_x_, thickness_);
  }
}

void Airfoil::ExtractCamberFromCoordinates(const std::vector<Scalar>& x_coords,
                                           const std::vector<Scalar>& y_coords) {
  // Detect format: Selig format has coordinates starting from TE upper,
  // going around LE to TE lower. X values decrease then increase.

  // Find the leading edge (minimum x)
  auto min_it = std::min_element(x_coords.begin(), x_coords.end());
  size_t le_index = std::distance(x_coords.begin(), min_it);

  // Split into upper and lower surfaces
  std::vector<Scalar> x_upper, y_upper, x_lower, y_lower;

  // Upper surface: from start to LE (Selig) or first half (Lednicer)
  for (size_t i = 0; i <= le_index; ++i) {
    x_upper.push_back(x_coords[i]);
    y_upper.push_back(y_coords[i]);
  }

  // Lower surface: from LE to end
  for (size_t i = le_index; i < x_coords.size(); ++i) {
    x_lower.push_back(x_coords[i]);
    y_lower.push_back(y_coords[i]);
  }

  // Reverse upper surface so it goes from LE to TE (increasing x)
  std::reverse(x_upper.begin(), x_upper.end());
  std::reverse(y_upper.begin(), y_upper.end());

  // Build splines for upper and lower surfaces
  CubicSpline upper_spline, lower_spline;
  if (x_upper.size() >= 2) {
    upper_spline.Build(x_upper, y_upper);
  }
  if (x_lower.size() >= 2) {
    lower_spline.Build(x_lower, y_lower);
  }

  // Generate uniform x distribution for camber line
  int num_points = 50;
  camber_x_.resize(num_points);
  camber_y_.resize(num_points);
  thickness_.resize(num_points);

  for (int i = 0; i < num_points; ++i) {
    Scalar x = static_cast<Scalar>(i) / static_cast<Scalar>(num_points - 1);
    camber_x_[i] = x;

    // Interpolate upper and lower surface
    Scalar y_up = upper_spline.IsValid() ? upper_spline.Evaluate(x) : 0.0;
    Scalar y_lo = lower_spline.IsValid() ? lower_spline.Evaluate(x) : 0.0;

    // Camber = mean line
    camber_y_[i] = 0.5 * (y_up + y_lo);

    // Thickness = distance between surfaces
    thickness_[i] = y_up - y_lo;
  }
}

}  // namespace vlm
