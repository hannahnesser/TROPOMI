'''
This is a script that takes as input 
lat_min, lat_max, lon_min, lon_max, name
and produces a netcdf from each of the
monthly regridded TROPOMI files.
'''

import pandas as pd
import numpy as np
from os import listdir, rename, mkdir
from os.path import join
import sys
sys.path.append('./python/')
import plot_oversampling as plot

def read_oversampling_output(file, regions, 
                             widths=[6,6,12,12,15,6], 
                             names=['row', 'col', 'lat', 'lon', 'xch4', 'cnt'], 
                             usecols=['lat', 'lon', 'xch4', 'cnt'],
                             chunksize=1e6):
    dfs = {r : pd.DataFrame(columns=usecols) for r in regions}
    for chunk in pd.read_fwf(file, widths=widths, 
                             header=None, names=names, usecols=usecols, 
                             chunksize=chunksize):
        for r_name, r_df in dfs.items():
            tmp = subset(chunk, regions[r_name])
            if tmp.shape[0] > 0:
                dfs[r_name] = dfs[r_name].append(tmp)
    return dfs

def subset(chunk, lat_lon_bnds):
    return chunk[(chunk['lat'] >= lat_lon_bnds[0]) &
                 (chunk['lat'] <= lat_lon_bnds[1]) &
                 (chunk['lon'] >= lat_lon_bnds[2]) &
                 (chunk['lon'] <= lat_lon_bnds[3])]

if __name__ == '__main__':
    data_dir = sys.argv[1]
    files = listdir(data_dir)
    files.sort()

    save_dir = sys.argv[2]
    
    regions = {}
    for r in sys.argv[3:]:
        regions[r.split(',')[0]] = np.array(r.split(',')[1:]).astype(float)    

    for f in files:
        print('')
        print('Processing ', f)
       
        # Check if the file name is incorrect
        if f.split('.')[-1] == 'csv_oversampled':
            print('Renaming input file.')
            fname = f.split('.')[0] + '_oversampled.csv'
            rename(join(data_dir, f), join(data_dir, fname))
            print('File renamed.')
       # If the file name ends in csv, we can proceed
        elif f.split('.')[-1] == 'csv':
            fname = f
        else:
            sys.exit('File name incorrect.')

        # Now we read in the file as a csv.
        print('Reading in data in chunks.')
        d = read_oversampling_output(join(data_dir, fname), regions)
        print('Data read.')
        print('Data subsetted.')

        # Now we make a directory to save our output.
        # First we check if the directory exits
        for r in regions.keys():
            if d[r].shape[0] == 0:
                print('No data remaining after subsetting.')
            else:
                if r not in listdir(save_dir):
                    print('Making save directory.')
                    mkdir(join(save_dir, r))
                # Ultimately, we want the file as a csv and netcdf. 
                # This will be easier to manipulate.
                save_name_date = f[:6] + '_' + r

                # csv
                d[r].to_csv(join(save_dir, r, save_name_date + '.csv'), index=False)
