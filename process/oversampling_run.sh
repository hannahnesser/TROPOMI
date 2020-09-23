#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 15000
#SBATCH -t 0-10:00
#SBATCH --mail-type=END
#SBATCH -J over_samp
#SBATCH -o slurm.%x.%j.out # STDOUT
#SBATCH -e slurm.%x.%j.err # STDERR

cd /n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/
cp ${HOME}/TROPOMI/fortran/*.* .
cp ${HOME}/TROPOMI/process/oversampling.sh .

./oversampling.sh
