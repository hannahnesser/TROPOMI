#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 600
#SBATCH -t 0-00:20
#SBATCH --mail-type=END
#SBATCH -J dams
#SBATCH -o slurm.%x.%j.out # STDOUT
#SBATCH -e slurm.%x.%j.err # STDERR

source activate invpy

# regions and lat/lon limits
# "region_name,latmin,latmax,lonmin,lonmax,lat,lon"
ITAIPU="itaipu,-27.9,-22.9,-57.1,-52.1,-25.40804,-54.58889"
TUCURUI="tucurui,-6.3,-1.3,-52.1,-47.1,-3.8323,-49.64926"
ROBERTBOURASSA="robertbourassa,51.3,56.3,-80.0,-75.0,53.77944,-77.54833"
ITUMBIARA="itumbiara,-20.9,-15.9,-51.6,-46.6,-18.40698,-49.09814"
PAMPULHA="pampulha,-22.3,-17.3,-46.5,-41.5,-19.846186,-43.966581"
SEGREDO="segredo,-28.3,-23.3,-54.6,-49.6,-25.793056,-52.113056"
HARTWELL="hartwell,31.9,36.9,-85.3,-80.3,34.357778,-82.821389"
FURNAS="furnas,-23.2,-18.2,-48.8,-43.8,-20.669722,-46.318056"
LAFORGE="laforge,51.7,56.7,-75.1,-70.1,54.166667,-72.616667"
BARRABONITA="barrabonita,-25.0,-20.0,-51.0,-46.0,-22.52029,-48.53565"
FUNIL="funil,-25.0,-20.0,-47.1,-42.1,-22.5287,-44.568"
BALBINA="balbina,-4.4,0.6,-62.0,-57.0,-1.917222,-59.473611"
WATTSBAR="wattsbar,33.1,38.1,-87.3,-82.3,35.621303,-84.781597"
EASTMAIN="eastmain,49.7,54.7,-78.4,-73.4,52.181667,-75.873889"
DOUGLAS="douglas,33.5,38.5,-86.0,-81.0,35.961111,-83.538889"
KARIBA="kariba,-19.0,-14.0,26.3,31.3,-16.522222,28.761667"

ALL=($ITAIPU $TUCURUI $ROBERTBOURASSA $ITUMBIARA $PAMPULHA $SEGREDO $HARTWELL $FURNAS $LAFORGE $BARRA $FUNIL $BALBINA $WATTSBAR $EASTMAIN $DOUGLAS $KARIBA)

# Code directory3
PYDIR=${HOME}/TROPOMI/python
INDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs_14_14/base/"
OUTDIR="/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs_14_14/dams/"

#python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
#python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
python ${PYDIR}/plot_oversampling_dams.py $OUTDIR $region
done
