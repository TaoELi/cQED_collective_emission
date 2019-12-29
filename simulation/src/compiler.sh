#!/bin/bash

# For QED solver, we complie and run it
g++ qed_solver.cpp -o run_qed -std=c++11 -I ~/local/armadillo-9.600.5/include -L ~/local/armadillo-9.600.5/lib64/ -larmadillo -O3 
wait
./run_qed
wait
# For MMST solver, we compile it with mpic++
mpic++ mmst_solver.cpp -o run_mmst -std=c++11 -I ~/local/armadillo-9.600.5/include -L ~/local/armadillo-9.600.5/lib64/ -larmadillo -O3 
wait
# On my personl computer
#mpirun -np 4 ./run_mmst
# On Nersc supercomputer
sbatch perform_nersc.sh
wait
