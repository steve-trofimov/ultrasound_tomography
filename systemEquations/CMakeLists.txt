add_library(systemEquations INTERFACE)

target_sources(systemEquations INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/systemEquations.h
				  ${CMAKE_CURRENT_SOURCE_DIR}/systemEquations.cpp)

target_include_directories(systemEquations INTERFACE
    "${PROJECT_SOURCE_DIR}/systemEquations/"
)

target_link_libraries(systemEquations INTERFACE Eigen3::Eigen ${CERES_LIBRARIES})

