#!/bin/bash

#SBATCH --job-name=normal2000
#SBATCH --nodes=16
#SBATCH --tasks-per-node=1
#SBATCH --mem=2GB 
#SBATCH --output=normal2000.out
#SBATCH --time=02:00:00 


module purge 
module load  intel/19.1.2 openmpi/intel/4.0.5

mpiexec ./partition normal2000.txt 2000
