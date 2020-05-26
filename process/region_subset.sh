#!/bin/bash

#SBATCH -n 15
#SBATCH -N 1
#SBATCH -p huce_amd
#SBATCH --mem 60000
#SBATCH -t 0-04:00
#SBATCH --mail-type=END

source activate invpy

# regions and lat/lon limits
# "region_name,latmin,latmax,lonmin,lonmax"
# Mexico coordinates:
# 14N to 33N Latitude, 86W to 118W Longitude
MEX="mexico,14,33,-118,-86"
ALL=($MEX)

# Code directory
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/regions/"

#python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done
