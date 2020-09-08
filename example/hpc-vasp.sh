#!/bin/bash
#SBATCH --job-name=jobname

#SBATCH --partition=etna
#SBATCH --account=nano
#SBATCH --qos=normal

#SBATCH --nodes=2
#SBATCH --time=01:00:00
#SBATCH --error=log-%A.err
#SBATCH --output=log-%A.job


## Run command
module purge
module load vasp_intelmpi

mpirun vasp_std > log
