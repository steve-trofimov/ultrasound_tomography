add_library(detector INTERFACE)

target_sources(detector INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/detector.h
				  ${CMAKE_CURRENT_SOURCE_DIR}/detector.cpp)

target_include_directories(detector INTERFACE
	"${PROJECT_SOURCE_DIR}/object_detection/"
)

