#!/bin/bash

#SBATCH -n 4
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 40000
#SBATCH -t 0-03:00
#SBATCH --mail-type=END
#SBATCH -J os_csv
#SBATCH -o slurm.%x.%j.out # STDOUT
#SBATCH -e slurm.%x.%j.err # STDERR

# download_subdir="/n/seasasfs02/hnesser/TROPOMI/downloads_14_14"
ROOTDIR="/n/holyscratch01/jacob_lab/hnesser/TROPOMI"
DATADIR="${ROOTDIR}/downloads_new/processed"
SAVEDIR="${ROOTDIR}/oversampling/input_csvs"
CODEDIR="${HOME}/TROPOMI"

# either "Lorente2020" or "operational"
data_version="Lorente2020"
MINDATE="20180101" # inclusive

# Make the savedir, if necessary
if [[ ! -d $SAVEDIR ]]; then
  mkdir -p $SAVEDIR
fi

# move to holyscratch
cd $ROOTDIR

# Actiate python environment
source activate troppy

# Run the correct script
if [[ data_version -eq "operational" ]]; then
  cp ${CODEDIR}/python/create_oversampling_csv_oper.py .
  python create_oversampling_csv_oper.py $DATADIR $SAVEDIR $MINDATE
elif [[ data_version -eq "Lorente2020" ]]; then
  cp ${CODEDIR}/python/create_oversampling_csv.py .
  python create_oversampling_csv.py $DATADIR $SAVEDIR $MINDATE
fi
