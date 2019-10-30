#!/bin/bash

cd build
mpirun -np 3 ./field_distance_matrix --N=1000 --tmin=2000 --tmax=10000 --tstep=2000 --nx=8 --ny=8 --nz=8 --inputfolder_relative=Liggghts/N_60000_U_1.333/ --inputfilenamepart=dump --outputfolder_relative=Liggghts/N_60000_U_1.333_1_5_2000_8_8_8/
cd ..
