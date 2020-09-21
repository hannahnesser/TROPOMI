#!/bin/bash
#
#SBATCH -p huce_intel
#SBATCH -n 2      # cores requested
#SBATCH -N 1      # nodes requested
#SBATCH --mem=2000  # memory in Mb
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 0-05:00  # time requested in hour:minute:second
#SBATCH --mail-type=END

# Use modules to set the software environment

#module load python/2.7.14-fasrc01
source activate invpy

export QT_QPA_PLATFORM='offscreen'

python ../python/Compare_TROPOMItoGC.py
