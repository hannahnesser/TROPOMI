#!/bin/bash

#list=(seq -f "%02g" 5 7)
list=("01" "06" "09" "11")
for i in ${list[@]}
do
sbatch ./download_data.sh $i
done
