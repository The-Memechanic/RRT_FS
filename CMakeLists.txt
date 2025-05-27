cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(RRT VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the SFML package
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)

# Add executables
add_executable(rrt rrt.cpp)

# Link the SFML libraries to both executables
target_link_libraries(rrt sfml-graphics sfml-window sfml-system)
