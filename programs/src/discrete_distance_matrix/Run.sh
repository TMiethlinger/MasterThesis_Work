#!/bin/bash

cd build
mpirun -np 3 ./discrete_distance_matrix --dim=6 --N=1000 --N_match=100 --l_inv=1 --v_inv=0 --tmin=2000 --tmax=10000 --tstep=2000 --inputfolder_relative=Liggghts/N_50000_U_1.333/ --inputfilenamepart=dump --outputfolder_relative=Liggghts/N_1000_U_1.333_1_5_2000_100/ --output_mode=1
cd ..
