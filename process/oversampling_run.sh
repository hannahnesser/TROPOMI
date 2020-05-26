#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 20000
#SBATCH -t 0-01:30
#SBATCH --mail-type=END

./create_oversampling_csv_run.sh
