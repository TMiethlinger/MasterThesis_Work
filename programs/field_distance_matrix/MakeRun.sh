#!/bin/bash



cd build
make
mpirun -np 3 ./field_distance_matrix --N=1000 --imin=2000 --imax=10000 --istep=2000 --nx=8 --ny=8 --nz=8 --inputfolder_relative=Liggghts/N_60000_U_1.333/ --inputfilenamepart=dump --outputfolder_relative=Liggghts/N_60000_U_1.333_1_5_2000_8_8_8/
