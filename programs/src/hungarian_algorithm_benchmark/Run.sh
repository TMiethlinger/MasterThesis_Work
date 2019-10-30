#!/bin/bash

cd build

Nlist=(100 200 400 600 800 1000 2000 4000 6000 8000 10000 20000 40000)

filenamepart="../../../../results/hungarian_algorithm_benchmark/"

filename="${filenamepart}timing_adjmat.txt"
for i in {1..3}
do
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        Nm=$N
        ./hungarian_algorithm_benchmark --N=$N --N_match=${Nm} --dim=6 >> ${filename}
    done
done

filename="${filenamepart}timing_adjlst_10.txt"
for i in {1..3}
do
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        Nm=$((N / 10))
        ./hungarian_algorithm_benchmark --N=$N --N_match=${Nm} --dim=6 >> ${filename}
    done
done

filename="${filenamepart}timing_adjlst_5.txt"
for i in {1..3}
do
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        Nm=$((N / 5))
        ./hungarian_algorithm_benchmark --N=$N --N_match=${Nm} --dim=6 >> ${filename}
    done
done

filename="${filenamepart}timing_adjlst_4.txt"
for i in {1..3}
do
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        Nm=$((N / 4))
        ./hungarian_algorithm_benchmark --N=$N --N_match=${Nm} --dim=6 >> ${filename}
    done
done

filename="${filenamepart}timing_adjlst_2.txt"
for i in {1..3}
do
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        Nm=$((N / 2))
        ./hungarian_algorithm_benchmark --N=$N --N_match=${Nm} --dim=6 >> ${filename}
    done
done
