#!/bin/bash

#SBATCH -n 15
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 60000
#SBATCH -t 0-08:00
#SBATCH --mail-type=END

source activate invpy

# regions and lat/lon limits
# "region_name,latmin,latmax,lonmin,lonmax"
VMA="vacamuerta,-40,-37,-70.5,-68.5"
SLG="sulige,37,40,107,110"
OBG="obagi,4.2,6.3,5.6,7.7"
KRK="kirkuk,34.8,36.9,43.0,45.1"
RLN="raslaffan,25,27,51,52"
#ALL=($VMA $SLG $OBG $RLN)
ALL=($KRK)

# Code directory
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/ong/"

python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
#python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done
