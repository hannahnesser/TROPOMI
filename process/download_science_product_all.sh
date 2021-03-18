#!/bin/bash

download_dir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new"
cp ./download_science_product.sh ${download_dir}
cd ${download_dir}

#for i in `seq -f "%02g" 0 15`


list=("07" "11")
for i in ${list[@]}
do
    sbatch download_science_product.sh $i
    sleep 30s
done
