#!/bin/bash

#PBS -o output/pbs.out
#PBS -j oe
#PBS -l nodes=10:ppn=8
#PBS -l walltime=500:00:00

module load gcc

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

###mpiexec ./prog -N 60000 -M 1000 -B 1000 -i 1 -j 200 -d 10000 > output/std_cout.txt
mpiexec ./fieldrecurrences -N 60000 -i 1 -j 12000 -d 1000

