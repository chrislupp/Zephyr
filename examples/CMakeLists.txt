# Examples CMakeLists.txt

# Basic analysis example
add_executable(basic_analysis basic_analysis.cpp)
target_link_libraries(basic_analysis PRIVATE vlm)

# Trim flight example
add_executable(trim_flight trim_flight.cpp)
target_link_libraries(trim_flight PRIVATE vlm)

# Stability analysis example
add_executable(stability_analysis stability_analysis.cpp)
target_link_libraries(stability_analysis PRIVATE vlm)
