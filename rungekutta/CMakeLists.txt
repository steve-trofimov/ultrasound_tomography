add_library(rungekutta INTERFACE)

target_sources(rungekutta INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/rungekutta.h)

target_include_directories(rungekutta INTERFACE
    "${PROJECT_SOURCE_DIR}/rungekutta/"
)

target_link_libraries(rungekutta INTERFACE Eigen3::Eigen)
