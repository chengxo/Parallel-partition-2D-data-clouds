#!/bin/bash

#SBATCH --job-name=normal8000
#SBATCH --nodes=64
#SBATCH --tasks-per-node=1
#SBATCH --mem=2GB 
#SBATCH --output=normal8000.out
#SBATCH --time=02:00:00 


module purge 
module load  intel/19.1.2 openmpi/intel/4.0.5

mpiexec ./partition normal8000.txt 8000
