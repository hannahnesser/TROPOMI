#!/bin/bash

#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 1000
#SBATCH -t 0-01:00:00

source activate invpy
python ~/TROPOMI/python/process_TROPOMI.py 3000 4080
