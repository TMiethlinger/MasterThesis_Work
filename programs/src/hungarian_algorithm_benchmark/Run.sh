#!/bin/bash

cd build

Nlist=(10 20 40 70 100 200 400 700 1000 2000 4000 7000 10000 20000 40000)

filenamepart="../../../../results/hungarian_algorithm_benchmark/"

for i in {1..5}
do
    filename="${filenamepart}timing_full_$i.txt"
    rm -rf ${filename}
    for N in ${Nlist[@]}; do
        ./hungarian_algorithm_benchmark --N=$N --dim=6 >> ${filename}
    done
done
