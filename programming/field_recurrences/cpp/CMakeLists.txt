#mkdir build
#cd build
#cmake ..
#make
cmake_minimum_required(VERSION 3.0)
project(field_recurrences CXX)

#BOOST
find_package(Boost 1.6 COMPONENTS program_options)

#MPI
find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_FLAGS  "-std=c++17 -Wall -O3")
SET(CMAKE_EXE_LINKER_FLAGS  "-lboost_program_options")

link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

add_executable(field_recurrences gemm_mpi_v1.cpp)
target_link_libraries(field_recurrences PUBLIC ${MPI_CXX_LIBRARIES} PUBLIC ${Boost_LIBRARIES})
