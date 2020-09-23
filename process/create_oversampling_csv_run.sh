#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 40000
#SBATCH -t 0-03:00
#SBATCH --mail-type=END

# download_subdir="/n/seasasfs02/hnesser/TROPOMI/downloads_14_14"
DATADIR="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/processed"
SAVEDIR="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/input_csvs"

# either "Lorente2020" or "operational"
data_version="Lorente2020"
MINDATE="20180101" # inclusive

# Make the savedir, if necessary
if [[ ! -d $output_dir ]]; then
  mkdir -p $output_dir
fi

# Actiate python environment
source activate troppy

# Run the correct script
if [[ data_version -eq "operational" ]]; then
  python create_oversampling_csv_oper.py $DATADIR $SAVEDIR $MINDATE
elif [[ data_version -eq "Lorente2020" ]]; then
  python create_oversampling_csv.py $DATADIR $SAVEDIR $MINDATE
fi
