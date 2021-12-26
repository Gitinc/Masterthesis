cmake_minimum_required(VERSION 3.21)
project(Masterarbeit)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Masterarbeit
        include/helper.h
        include/LorentzAttractor.h
        include/MackeyGlass.h
        include/Mask.h
        include/Narma10.h
        src/helper.cpp
        src/LorentzAttractor.cpp
        src/MackeyGlass.cpp
        src/Mask.cpp
        src/Narma10.cpp
        src/TestSL.cpp
        main.cpp)
target_link_libraries(Masterarbeit ${ARMADILLO_LIBRARIES})