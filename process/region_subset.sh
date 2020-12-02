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
ITAIPU="ITAIPU,-26.4,-24.4,-55.6,-53.6"
TUCURUI="TUCURUI,-4.9,-2.9,-50.6,-48.6"
ROBERTBOURASSA="ROBERTBOURASSA,52.7,54.7,-78.0,-76.0"
ITUMBIARA="ITUMBIARA,-19.4,-17.4,-50.1,-48.1"
PAMPULHA="PAMPULHA,-20.8,-18.8,-45.0,-43.0"
SEGREDO="SEGREDO,-26.8,-24.8,-53.1,-51.1"
HARTWELL="HARTWELL,33.5,35.5,-83.9,-81.9"
FURNAS="FURNAS,-21.7,-19.7,-47.3,-45.3"
LAFORGE="LAFORGE,53.2,55.2,-73.6,-71.6"
BARRA="BARRA,BONITA-23.5,-21.5,-49.6,-47.6"
FUNIL="FUNIL,-23.5,-21.5,-45.6,-43.6"
BALBINA="BALBINA,-2.9,-0.9,-60.5,-58.5"
WATTSBAR="WATTSBAR,34.6,36.6,-85.8,-83.8"
EASTMAIN="EASTMAIN,50.5,52.5,-75.0,-73.0"
DOUGLAS="DOUGLAS,35.0,37.0,-84.5,-82.5"
KARIBA="KARIBA,-18.0,-16.0,27.0,29.0"
ALL=($ITAIPU $TUCURUI $ROBERTBOURASSA $ITUMBIARA $PAMPULHA $SEGREDO $HARTWELL $FURNAS $LAFORGE $BARRA $FUNIL $BALBINA $WATTSBAR $EASTMAIN $DOUGLAS $KARIBA)

# Code directory
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/dams/"

python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
#python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done
