#!/bin/bash -l

#SBATCH -q regular
#SBATCH --nodes=128
#SBATCH --tasks-per-node=64
#SBATCH -t 1:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C knl
#SBATCH -A m3138
 

srun --cpu-bind=cores ./run_mmst

