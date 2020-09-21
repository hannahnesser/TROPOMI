#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 1000
#SBATCH -t 0-06:00

python ${1}/process_science_product.py ${2} ${3} ${4} "${@:5}"
