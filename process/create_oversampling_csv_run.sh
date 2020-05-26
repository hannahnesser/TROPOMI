#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 40000
#SBATCH -t 0-03:00
#SBATCH --mail-type=END

download_subdir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads"
output_dir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling_input_csvs"
MINDATE="20180401" # inclusive

# Make the savedir, if necessary


source activate invpy

python create_oversampling_csv.py $DATADIR $SAVEDIR $MINDATE
