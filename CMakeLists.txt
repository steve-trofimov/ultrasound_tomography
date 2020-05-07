cmake_minimum_required(VERSION 3.10)

project(timeOfFlight)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Ceres REQUIRED)

add_subdirectory(object_detection)
add_subdirectory(oneToF)
add_subdirectory(read_data)
add_subdirectory(systemEquations)
add_subdirectory(rungekutta)
add_subdirectory(app)
