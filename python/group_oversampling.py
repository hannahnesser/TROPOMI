from os import listdir
from os.path import join
#import xarray as xr
import pandas as pd
import sys

data_dir = sys.argv[1]
region = sys.argv[2]
data_dir = join(data_dir, region)

files = listdir(data_dir)
files = [f for f in files if f[-3:] == 'csv']
files.sort(reverse=True)

yr_mo = [f[:6] for f in files]
yrs = [int(f[:4]) for f in files]
yrs = list(set(yrs))
yrs.sort(reverse=True)

seasons = {}
for y in yrs:
    seasons[str(y) + 'MAM'] = [str(y) + m for m in ['03', '04', '05']]
    seasons[str(y) + 'JJA'] = [str(y) + m for m in ['06', '07', '08']]
    seasons[str(y) + 'SON'] = [str(y) + m for m in ['09', '10', '11']]
    if (y-1) in yrs:
        seasons[str(y)  + 'DJF'] = [str(y-1) + '12', str(y) + '01', str(y) + '02']

for s, mlist in seasons.items():
    if set(mlist).issubset(yr_mo):
        s_files = [f for f in files if f[:6] in mlist]
        s_csvs = []
        for sf in s_files:
            d = pd.read_csv(join(data_dir, sf))
            d['xch4'] *= d['cnt']
            d = d.set_index(['lat', 'lon'])
            s_csvs.append(d)
        t = pd.concat(s_csvs, axis=1, join='outer')
        t['cnt_tot'] = t['cnt'].sum(axis=1)
        t['xch4_tot'] = t['xch4'].sum(axis=1)/t['cnt_tot']
        t = t.drop(columns=['xch4', 'cnt'])
        t = t.rename(columns={'xch4_tot' : 'xch4', 'cnt_tot' : 'cnt'})
        t = t.reset_index()
        t.to_csv(join(data_dir, s + '_' + region + '.csv'), index=False)
        
