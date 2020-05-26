#!/bin/bash

#SBATCH --mail-type=END

## Define our variables

# Download directories (i.e. where to download the raw data)
download_dir="/n/holyscratch01/jacob_lab/hnesser/TROPOMI"
download_subdir="downloads"

# Script directories (i.e. where to find the original scripts)
script_dir="/n/home04/hnesser/TROPOMI"

# Last download information (i.e. information on the last set of
# downloads) 
last_download_dir="/n/seasasfs02/hnesser/TROPOMI"
last_download_name="lastdownload"
last_download=${last_download_dir}/${last_download_name}

## Then download the new data

# First, download all new TROPOMI data
sbatch download_run.sh ${download_dir} ${download_subdir} ${script_dir} ${last_download_dir} ${last_download_name}

# Second, create the oversampling csvs
