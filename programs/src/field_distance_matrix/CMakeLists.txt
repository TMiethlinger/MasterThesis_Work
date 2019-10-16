#mkdir build
#cd build
#cmake -DBoost_INCLUDE_DIR="/home/k3501/k354524/apps/boost_1_71_0/" -DBoost_LIBRARY_DIR="/home/k3501/k354524/apps/boost_1_71_0/stage/lib/" ..
#make
#./helloworld --n=...
cmake_minimum_required(VERSION 3.0)
project(field_distance_matrix CXX)

#Boost
FIND_PACKAGE(Boost 1.6 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#MPI
find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CXX_FLAGS  "-std=c++17 -Wall -O3")
SET(CMAKE_EXE_LINKER_FLAGS  "-lboost_program_options")

add_executable(field_distance_matrix field_distance_matrix.cpp ../../includes/general_util.hpp ../../includes/io_util.hpp)
target_link_libraries(field_distance_matrix PUBLIC ${MPI_CXX_LIBRARIES} PUBLIC ${Boost_LIBRARIES})
