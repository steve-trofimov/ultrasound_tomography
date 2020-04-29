file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(app ${SRCS} ../object_detection/detector.cpp ../object_detection/detector.h argparse.h)

find_package(Boost COMPONENTS program_options REQUIRED)

TARGET_LINK_LIBRARIES(app oneToF read_data detector ${Boost_LIBRARIES})
