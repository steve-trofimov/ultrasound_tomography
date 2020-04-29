add_library(oneToF INTERFACE)

target_sources(oneToF INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/oneToF.h)

target_include_directories(oneToF INTERFACE
	"${PROJECT_SOURCE_DIR}/oneToF/"
)
