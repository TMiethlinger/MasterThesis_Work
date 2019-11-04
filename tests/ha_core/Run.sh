#!/bin/bash

cd build

N=2000
./ha_core_test --N=$N --N_match=450 --dim=6

