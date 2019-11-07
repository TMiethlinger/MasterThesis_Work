#!/bin/bash

cd build

N=1000
./ha_core_test --N=$N --N_match=300 --dim=6
