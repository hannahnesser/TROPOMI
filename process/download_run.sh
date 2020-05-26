#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 200
#SBATCH -t 0-20:00
#SBATCH --mail-type=END
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

download_dir=$1
download_subdir=$2
script_dir=$3
last_download_dir=$4
last_download_name=$5

./download.sh ${download_dir} ${download_subdir} ${script_dir} ${last_download_dir} ${last_download_name}
