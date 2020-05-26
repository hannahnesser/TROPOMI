#!/bin/bash

#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 1000
#SBATCH -t 0-0:60
#SBATCH --mail-type=END

source activate invpy
python process_TROPOMI.py
