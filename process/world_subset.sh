#!/bin/bash

#SBATCH -n 3
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 10000
#SBATCH -t 0-10:00
#SBATCH --mail-type=END
#SBATCH -J world
#SBATCH -o slurm.%x.%j.out # STDOUT
#SBATCH -e slurm.%x.%j.err # STDERR

source activate invpy
export OMP_NUM_THREADS=4

# regions and lat/lon limits
# "region_name,latmin,latmax,lonmin,lonmax"
# NAM="northamerica,15,90,-180,-18"
# SAM="southamerica,-60,15,-105,-30"
# AFR="africa,-35,36,-18,60"
# EUR="eurasia,36,90,-18,180"
# PAC="pacific,-53,36,60,180"
# ALL=($NAM $SAM $AFR $EUR $PAC)
#ALL=($NAM $EUR $PAC)

NAM="northamerica,10.375,59.375,-129.21875,-60.78125"
ALL=($NAM)

# Code directory
# PYDIR=${HOME}/TROPOMI/python
# INDIR="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/output_csvs/"
# OUTDIR="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/oversampling/output_csvs/world/"

BASE="${HOME}/Documents/Harvard/Research/TROPOMI_Inversion"
PYDIR="${BASE}/python"
INDIR="${BASE}/observations/oversampling"
OUTDIR="${BASE}/observations/oversampling"

if [[ ! -d $OUTDIR ]]; then
  mkdir -p $OUTDIR
fi

# Only needs to be run once.
#python ${PYDIR}/subset_oversampling.py $INDIR $OUTDIR "${ALL[@]}"

for region in "${ALL[@]}"
do
echo $region
#python ${PYDIR}/group_oversampling.py $OUTDIR "${region%%,*}"
python ${PYDIR}/plot_oversampling.py $OUTDIR $region
done

#cp ${INDIR}*.csv* /n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs_14_14/base/
#cp -r ${OUTDIR} /n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs_14_14/world
