#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 20000
#SBATCH -t 0-01:30
#SBATCH --mail-type=END
#SBATCH -J over_samp
#SBATCH -o slurm.%x.%j.out # STDOUT
#SBATCH -e slurm.%x.%j.err # STDERR

./oversampling.sh
