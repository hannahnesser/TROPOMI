#!/bin/bash

#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 1000
#SBATCH -t 0-01:00

cd raw
wget ftp://ftp.sron.nl/open-access-data-2/TROPOMI/tropomi/ch4/14_14_Lorente_et_al_2020_AMTD/s5p_l2_ch4_0014_$1*
