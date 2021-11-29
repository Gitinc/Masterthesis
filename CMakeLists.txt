cmake_minimum_required(VERSION 3.20)
project(untitled)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD)

add_executable(untitled main.cpp helper.cpp helper.h MackeyGlass.cpp TestSL.cpp MackeyGlass.h Mask.cpp Mask.h)
target_link_libraries(untitled ${ARMADILLO_LIBRARIES})

