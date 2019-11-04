#!/bin/bash

#PBS -o pbs/pbs_output.txt
#PBS -e pbs/pbs_error.txt
#PBS -l nodes=8:ppn=8
#PBS -l walltime=24:00:00

module load gcc
module load mvapich2

cd $PBS_O_WORKDIR
cd build/

simset="Liggghts"
N="50000"
param="U_1.333"
dataset="N_${N}_${param}"
tstep=2000
tmin=${tstep}
tmax=12000000
imin=$((tmin / tstep))
imax=$((tmax / tstep))

n=12
mpiexec ./field_distance_matrix --N=${N} --tmin=${tmin} --tmax=${tmax} --tstep=${tstep} --nx=${n} --ny=${n} --nz=${n} --inputfolder_relative="${simset}/${dataset}/" --inputfilenamepart=dump --outputfolder_relative="${simset}/${dataset}_${imin}_${imax}_${tstep}_${n}_${n}_${n}/"
