#!/bin/bash

# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=SE_101TLS_FDTD_dx_13

mkdir $NEW_PATH
cp $ORIGIN_PATH/* $NEW_PATH

# Under new path, we modify some basic parameters
NTRAJS=12800

cd $NEW_PATH

sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 128/NTRAJ_MULTITRAJ_Eh = $NTRAJS/g" parameters.hpp
sed -i".bak" "s/NMOLECULES = 1/NMOLECULES = 101/g" parameters.hpp
sed -i".bak" "s/LATTICE_SPACING_GRIDS = 25/LATTICE_SPACING_GRIDS = 13/g" parameters.hpp
#sed -i".bak" "s/\/\/ #define RUN_MODES_INSTEAD_OF_FDTD/#define RUN_MODES_INSTEAD_OF_FDTD/g" parameters.hpp


./compiler.sh

cd ..
