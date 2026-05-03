# Examples CMakeLists.txt

# Basic analysis example
add_executable(basic_analysis basic_analysis/basic_analysis.cpp)
target_link_libraries(basic_analysis PRIVATE vlm)

# Trim flight example
add_executable(trim_flight trim_flight/trim_flight.cpp)
target_link_libraries(trim_flight PRIVATE vlm)

# Stability analysis example
add_executable(stability_analysis stability_analysis/stability_analysis.cpp)
target_link_libraries(stability_analysis PRIVATE vlm)

# VTK visualization example (builds even without VTK - checks availability at runtime)
add_executable(vtk_visualization vtk_visualization/vtk_visualization.cpp)
target_link_libraries(vtk_visualization PRIVATE vlm)

# Analysis scaling benchmark
add_executable(benchmark_analysis benchmark_analysis_scaling/benchmark_analysis.cpp)
target_link_libraries(benchmark_analysis PRIVATE vlm)
