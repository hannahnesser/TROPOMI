#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 1000
#SBATCH -t 0-06:00

## -------------------------------------------------------------------------##
## Load the environment
## -------------------------------------------------------------------------##
echo "Activating python environment"

module load Anaconda3/5.0.1-fasrc01
source activate ~/python/miniconda/envs/TROPOMI_inversion

echo "Activated ${CONDA_PREFIX}"

## -------------------------------------------------------------------------##
## Run the script
## -------------------------------------------------------------------------##
python ${1}/process_science_product.py ${2} ${3} ${4} "${@:5}"
