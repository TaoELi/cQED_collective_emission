#!/bin/bash -l

#SBATCH -q premium
#SBATCH --nodes=128
#SBATCH --tasks-per-node=64
#SBATCH -t 2:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C knl
#SBATCH -A m3138

# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=SE_Dicke

mkdir $NEW_PATH

# Under new path, we modify some basic parameters
NTRAJS=12800

cd $NEW_PATH

for N in 1 3 5 7 9 11 13 15 17 19 21 25 35 45

do
    CURRENT_SUB_PATH=N_$N
    mkdir $CURRENT_SUB_PATH
    cp ../$ORIGIN_PATH/* $CURRENT_SUB_PATH
    cd $CURRENT_SUB_PATH
    sed -i".bak" "s/C.row(1).zeros()/C.row(1).ones()/g" mmst_solver.cpp
    sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 128/NTRAJ_MULTITRAJ_Eh = $NTRAJS/g" parameters.hpp
    sed -i".bak" "s/NMOLECULES = 1/NMOLECULES = $N/g" parameters.hpp
    sed -i".bak" "s/1000/5000/g" parameters.hpp
    sed -i".bak" "s/(NGRIDS - 1)/(NGRIDS - 1)\/ 5.0/g" parameters.hpp
    sed -i".bak" "s/2500/2500, 2500, 2500/g" parameters.hpp
    sed -i".bak" "s/{1.0}/{1.0, 1.0, 1.0}/g" parameters.hpp
    sed -i".bak" "s/LATTICE_SPACING_GRIDS = 25/LATTICE_SPACING_GRIDS = 0/g" parameters.hpp
    #sed -i".bak" "s/\/\/ #define RUN_MODES_INSTEAD_OF_FDTD/#define RUN_MODES_INSTEAD_OF_FDTD/g" parameters.hpp
    
    if [[ ! -e "traj_MultiEh.txt" ]]
    then
        pwd
        echo "Compiling"
        #mpic++ mmst_solver.cpp -o run_mmst -std=c++11 -I ~/local/armadillo-9.600.5/include -L ~/local/armadillo-9.600.5/lib64/ -larmadillo -O3
        wait
        echo "End Compiling"
        srun --cpu-bind=cores ./run_mmst
        wait
    fi

    cd ..
done

cd ..
