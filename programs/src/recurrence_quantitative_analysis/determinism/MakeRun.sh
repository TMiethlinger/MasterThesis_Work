#!/bin/bash

cd build
make
mpirun -np 3 ./determinism --n=5 --neps=20 --inputfolder_relative=field_distance_matrix/Liggghts/N_60000_U_1.333_1_5_2000_8_8_8/ --inputfilename=field_distance_matrix.txt --outputfolder_relative=Liggghts/N_60000_U_1.333_1_5_2000_8_8_8/ --outputfilename=determinism_20.txt
