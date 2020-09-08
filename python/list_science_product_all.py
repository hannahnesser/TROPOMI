#!/usr/bin/env python3

import xarray as xr
import pandas as pd
from datetime import date
from os import listdir
from os.path import join
import sys
import process_science_product as psp

# I will only process orbits beginning with 2832 (May 1, 2018)
raw_data_dir = str(sys.argv[1])
first_orbit = int(sys.argv[2])

# Create a list of files
files = listdir(raw_data_dir)
files = [f for f in files if (f[-2:] == 'nc') and (len(f) == 24)] 
files = [f for f in files if int(f.split('_')[-1].split('.')[0]) >= first_orbit]
files.sort()

print(*files, sep=' ')
