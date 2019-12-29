#!/bin/bash

# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=EM_FDTD_Full

mkdir $NEW_PATH
cp $ORIGIN_PATH/* $NEW_PATH

# Under new path, we modify some basic parameters
NTRAJS=1280000

cd $NEW_PATH

# Modifying parameters.hpp
sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 128/NTRAJ_MULTITRAJ_Eh = $NTRAJS/g" parameters.hpp
#sed -i".bak" "s/\/\/ #define RUN_MODES_INSTEAD_OF_FDTD/#define RUN_MODES_INSTEAD_OF_FDTD/g" parameters.hpp

# Modifying nersc submit script
sed -i".bak" "s/1:00:00/10:00:00/g" perform_nersc.sh

# Compile
./compiler.sh

cd ..
