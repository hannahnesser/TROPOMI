#!/bin/bash

download_dir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/raw"

list=("10")

#for i in `seq -f "%02g" 1 15`

for i in ${list[@]}
do
    sbatch ./download_science_product.sh $i $download_dir
    sleep 30s
done
