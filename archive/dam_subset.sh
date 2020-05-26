#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 4000
#SBATCH -t 0-04:00
#SBATCH --mail-type=END

source activate invpy

# regions and lat/lon limits
# "region_name,latmin,latmax,lonmin,lonmax"
ITA="itaipu,-26.38,-24.38,-55.55,-53.55"
TUC="tucurui,-4.92,-2.92,-50.6,-48.6"
ROB="robertbourassa,52.67,54.67,-78,-76"
ITU="itumbiara,-19.4,-17.4,-50.05,-48.05"
PAM="pampulha,-20.83,-18.83,-44.97,-42.97"
SEG="segredo,-26.78,-24.78,-53.1,-51.1"
HAR="hartwell,33.47,35.47,-83.85,-81.85"
FUR="furnas,-21.67,-19.67,-47.25,-45.25"
LAF="laforge,53.2,55.2,-73.58,-71.58"
BAR="barrabonita,-23.52,-21.52,-49.58,-47.58"
FUN="funil,-23.52,-21.52,-45.55,-43.55"
BAL="balbina,-2.92,-0.92,-60.47,-58.47"
WAT="wattsbar,34.62,36.62,-85.78,-83.78"
EAS="eastmain,50.5,52.5,-75,-73"
DOU="douglas,34.96,36.96,-84.54,-82.54"
KAR="kariba,-18,-16,27,29"
ALL=($ITA $TUC $ROB $ITU $PAM $SEG $HAR $FUR $LAF $BAR $FUN $BAL $WAT $EAS $DOU $KAR)

# Code directory
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs/dams/"

python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done
