#!/bin/bash

# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=GroundState_Stability

mkdir $NEW_PATH
cp $ORIGIN_PATH/* $NEW_PATH

# Under new path, we modify some basic parameters
NTRAJS=12800

cd $NEW_PATH

sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 128/NTRAJ_MULTITRAJ_Eh = $NTRAJS/g" parameters.hpp
sed -i".bak" "s/EXCITEDSTATEDISTRIBUTION\[\] = {1.0}/EXCITEDSTATEDISTRIBUTION\[\] = {0.0}/g" parameters.hpp
#sed -i".bak" "s/\/\/ #define RUN_MODES_INSTEAD_OF_FDTD/#define RUN_MODES_INSTEAD_OF_FDTD/g" parameters.hpp

# Do turn off ZPE 
sed -i".bak" "s/\/\/EhrenfestMulti model_SED/EhrenfestMulti model_SED/g" mmst_solver.cpp
sed -i".bak" "s/\/\/model_SED.run_parallel/model_SED.run_parallel/g" mmst_solver.cpp
sed -i".bak" "s/\/\/EhrenfestMulti model_SQC/EhrenfestMulti model_SQC/g" mmst_solver.cpp
sed -i".bak" "s/\/\/model_SQC.run_parallel/model_SQC.run_parallel/g" mmst_solver.cpp


./compiler.sh

cd ..
