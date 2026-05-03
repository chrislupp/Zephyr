# Project Pegasus
# Version 0.1.0
#
# Copyright 2024-2026, Christopher A. Lupp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Examples CMakeLists.txt

# Basic analysis example
add_executable(basic_analysis basic_analysis/basic_analysis.cpp)
target_link_libraries(basic_analysis PRIVATE pegasus)

# Trim flight example
add_executable(trim_flight trim_flight/trim_flight.cpp)
target_link_libraries(trim_flight PRIVATE pegasus)

# Stability analysis example
add_executable(stability_analysis stability_analysis/stability_analysis.cpp)
target_link_libraries(stability_analysis PRIVATE pegasus)

# VTK visualization example (builds even without VTK - checks availability at runtime)
add_executable(vtk_visualization vtk_visualization/vtk_visualization.cpp)
target_link_libraries(vtk_visualization PRIVATE pegasus)

# Analysis scaling benchmark
add_executable(benchmark_analysis benchmark_analysis_scaling/benchmark_analysis.cpp)
target_link_libraries(benchmark_analysis PRIVATE pegasus)
