#!/bin/bash

# Note the raw data dir and the processed data dir
base=/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new
python_dir=${HOME}/TROPOMI/python
raw=${base}/raw
processed=${base}/processed
now="$(date +'%Y%m%dT%H%M')"

# Note the number of the first orbit from which to start
first_orbit=2832

# Now cd into base
cp process_science_product.sh $base
cd $base

# activate environment
source activate troppy

# Get list of files
files=($(python ${python_dir}/list_science_product_all.py $raw $first_orbit))
num=${#files[@]}

# pseudo code from here on out
i=0
n=500
while [ $i -lt $num ]
do
    file_subset="${files[@]:$i:$n}"
    sbatch process_science_product.sh $python_dir $raw $processed $now $file_subset
    #python ${python_dir}/process_science_product.py $raw $processed $now $file_subset
    i=$(($i + $n))
done
