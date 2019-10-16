#!/bin/bash

rm -rf build
mkdir build
cd build
cmake -DBoost_INCLUDE_DIR="/home/k3501/k354524/apps/boost_1_71_0/" -DBoost_LIBRARY_DIR="/home/k3501/k354524/apps/boost_1_71_0/stage/lib/" ..
cd ..
