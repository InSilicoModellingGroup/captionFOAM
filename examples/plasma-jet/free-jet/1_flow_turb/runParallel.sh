#!/bin/bash
#SBATCH --job-name=1.0_350intense
#SBATCH --output=mpijob.%x_%j.out
#SBATCH --error=mpijob.%x_%j.err
#SBATCH --ntasks=240
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=COMPUTE
#SBATCH --time=4-00:00:00
 
mpirun -np ${SLURM_NTASKS} rhoReactingBuoyantFoam -parallel