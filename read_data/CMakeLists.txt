add_library(read_data INTERFACE)

target_sources(read_data INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/read_data.h)

target_include_directories(read_data INTERFACE
	"${PROJECT_SOURCE_DIR}/read_data/"
)
