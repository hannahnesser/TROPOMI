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
NAM="northamerica,15,90,-180,-18"
SAM="southamerica,-60,15,-105,-30"
AFR="africa,-35,36,-18,60"
EUR="eurasia,36,90,-18,180"
PAC="pacific,-53,36,60,180"
ALL=($NAM $SAM $AFR $EUR $PAC)
#ALL=($NAM $EUR $PAC)

# Code directory
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/world/"

python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done
