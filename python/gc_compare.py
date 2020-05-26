import pandas as pd
import numpy as np
from os import listdir, rename, mkdir, environ
from os.path import join

# Plotting
from matplotlib import rcParams, colorbar, colors, cm
import matplotlib.pyplot as plt

# Tell matplotlib not to look for an X-window
environ['QT_QPA_PLATFORM']='offscreen'

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
#rcParams['text.usetex'] = True

def color(k):
    c = cm.get_cmap('inferno', lut=10)
    return colors.to_hex(c(k))

data_dir='/n/scratchlfs/jacob_lab/hnesser/GC_12.4.0/Global_4x5/J_Emis_Clean/run_dirs/J_emis_0000'
yr = 2018
mo = 7
lat_bands = np.arange(-90, 91, 30)

files = [f for f in listdir(data_dir)
         if (f[:15] == 'sat_obs.tropomi')]
files.sort()

df = pd.DataFrame()
for f in files:
    df = df.append(pd.read_csv(join(data_dir, f), delim_whitespace=True, header=0))

df['LAT_BAND'] = pd.cut(df['LAT'], bins=lat_bands)
df['TROPO'] *= 1e9
df['model'] *= 1e9

fig, ax = plt.subplots(2, 3, figsize=(15, 10))
plt.subplots_adjust(hspace=0.3, wspace=0.3)
lat_bands = np.unique(df['LAT_BAND'])
i = 0
for axis in ax.flatten()[:-1]:
    d = df[df['LAT_BAND'] == lat_bands[i]]
    axis.set_title(r'%s' % str(lat_bands[i]), fontsize=20)
    axis.scatter(d['TROPO'], d['model'], 
                 s=5, c=color(5), alpha=0.3)
    axis.set_xlabel('TROPOMI')
    axis.set_ylabel('Model')
    axis.set_xlim(axis.get_ylim())
    axis.set_ylim(axis.get_ylim())
    axis.plot(axis.get_ylim(), axis.get_ylim(), 
              c=color(1), ls='--', lw=1)
    axis.set_aspect('equal')
    i += 1

ax[1,2].axis('off')

fig.savefig(join(data_dir, 'gc_trop_comparison.png'))
