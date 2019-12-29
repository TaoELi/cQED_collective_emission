#!/bin/bash

# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=SE_1TLS_FDTD_near_mirror

mkdir $NEW_PATH

# Under new path, we modify some basic parameters
NTRAJS=12800

cd $NEW_PATH

for dx in 0 #13 25 50

do
    CURRENT_SUB_PATH=dx_$dx
    mkdir $CURRENT_SUB_PATH
    cp ../$ORIGIN_PATH/* $CURRENT_SUB_PATH
    cd $CURRENT_SUB_PATH
    sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 128/NTRAJ_MULTITRAJ_Eh = $NTRAJS/g" parameters.hpp
    sed -i".bak" "s/2500/$dx/g" parameters.hpp
    #sed -i".bak" "s/\/\/ #define RUN_MODES_INSTEAD_OF_FDTD/#define RUN_MODES_INSTEAD_OF_FDTD/g" parameters.hpp
    ./compiler.sh
    cd ..
done

cd ..
