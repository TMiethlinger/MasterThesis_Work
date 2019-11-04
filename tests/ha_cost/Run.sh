#!/bin/bash

cd build

N=4
./ha_cost_test --N=$N --N_match=2 --dim=6
