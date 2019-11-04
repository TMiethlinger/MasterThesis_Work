cmake_minimum_required(VERSION 3.0)
project(ha_core_test CXX)

#Boost
FIND_PACKAGE(Boost 1.6 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#MPI
find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_FLAGS  "-std=c++17 -Wall -O3")
SET(CMAKE_EXE_LINKER_FLAGS  "-lboost_program_options")

add_executable(ha_core_test ha_core_test.cpp ../../programs/src/hungarian_algorithm/ha_core.hpp ../../programs/src/hungarian_algorithm/ha_cost.hpp ../../programs/src/hungarian_algorithm/ha_distance.hpp ../../programs/includes/io_util.hpp )
target_link_libraries(ha_core_test PUBLIC ${MPI_CXX_LIBRARIES} PUBLIC ${Boost_LIBRARIES})
