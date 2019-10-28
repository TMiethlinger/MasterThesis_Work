cmake_minimum_required(VERSION 3.0)
project(hungarian_algorithm_benchmark CXX)

#Boost
FIND_PACKAGE(Boost 1.6 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#MPI
find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_FLAGS  "-std=c++17 -Wall -O3")
SET(CMAKE_EXE_LINKER_FLAGS  "-lboost_program_options")

add_executable(hungarian_algorithm_benchmark hungarian_algorithm_benchmark.cpp ../hungarian_algorithm/ha_core.hpp ../hungarian_algorithm/ha_costobject.hpp ../hungarian_algorithm/ha_distance.hpp)
target_link_libraries(hungarian_algorithm_benchmark PUBLIC ${MPI_CXX_LIBRARIES} PUBLIC ${Boost_LIBRARIES})
