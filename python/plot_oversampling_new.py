from os import listdir, mkdir, getcwd
from os.path import join, dirname, realpath
import sys

import pandas as pd
import numpy as np
import xarray as xr

import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
sys.path.append('.')
import format_plots as fp
import config as config

# from fastkml import kml
# from lxml import etree
# from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
#                        AltitudeMode, Camera)

import cartopy.crs as ccrs
import cartopy

import shapefile as shp

colors = plt.cm.get_cmap('inferno', lut=8)

# Other font details
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = config.LABEL_FONTSIZE*config.SCALE
rcParams['text.usetex'] = True

def plot_TROPOMI(data, latlim, lonlim, res, figax=None, title='', genkml=False, vals='xch4', **plot_options):
    # Set default plot options]
    if 'vmin' not in plot_options:
        plot_options['vmin'] = 1700
    if 'vmax' not in plot_options:
        plot_options['vmax'] = 1900
    if 'cmap' not in plot_options:
        plot_options['cmap'] = 'inferno'

    # Define edges
    lat_steps = (latlim[1]-latlim[0])/res + 1
    lon_steps = (lonlim[1]-lonlim[0])/res + 1
    if (lat_steps != int(lat_steps)) or (lon_steps != int(lon_steps)):
        sys.exit('Bad interval.')

    lats = np.around(np.linspace(latlim[0], latlim[1], int(lat_steps)) + res/2,
                     len(str(res/2).split('.')[1]))
    lats_s = pd.DataFrame({'idx' : np.ones(len(lats)-1), 'lat' : lats[:-1]})
    lons = np.around(np.linspace(lonlim[0], lonlim[1], int(lon_steps)) + res/2,
                     len(str(res/2).split('.')[1]))
    lons_s = pd.DataFrame({'idx' : np.ones(len(lons)-1), 'lon' : lons[:-1]})
    df = pd.merge(lats_s, lons_s, on='idx').drop(columns='idx')
    data = pd.merge(df, data, on=['lat', 'lon'], how='left')

    d_p = data.pivot(index='lat', columns='lon', values=vals)
    lon, lat = np.meshgrid(lons-res/2, lats-res/2)

    if figax is None:
        fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection' : ccrs.PlateCarree()})
    else:
        fig, ax = figax

    c = ax.pcolormesh(lons, lats, d_p,
                      snap=True,
                      vmin=plot_options['vmin'], vmax=plot_options['vmax'],
                      cmap=plot_options['cmap'],
                      edgecolors=None)

    return fig, ax, c

if __name__ == '__main__':
    DATA_DIR = sys.argv[1]
    REGION = str(sys.argv[2].split(',')[0])
    DATA_DIR = join(DATA_DIR, REGION)

    region_dict = {'ITAIPU' : 'Itaipu',
                   'TUCURUI' : 'Tucurui',
                   'ROBERTBOURASSA' : 'Robert-Bourassa (La Grande-2)',
                   'ITUMBIARA' : 'Itumbiara',
                   'PAMPULHA' : 'Pampulha',
                   'SEGREDO' : 'Segredo',
                   'HARTWELL' : 'Hartwell',
                   'FURNAS' : 'Furnas',
                   'LAFORGE' : 'Laforge-1',
                   'BARRABONITA' : 'Barra Bonita ',
                   'FUNIL' : 'Funil',
                   'BALBINA' : 'Balbina',
                   'WATTSBAR' : 'Watts Bar',
                   'EASTMAIN' : 'Eastmain 1',
                   'DOUGLAS' : 'Douglas',
                   'KARIBA' : 'Kariba'}

    latlim = np.array(sys.argv[2].split(',')[1:3])
    lonlim = np.array(sys.argv[2].split(',')[3:])
    res = 0.01
    vmin = 1800
    vmax = 1900
    count_min = 1
    plot_options = {'vmin' : vmin, 'vmax' : vmax, 'cmap' : 'inferno'}

    files = listdir(DATA_DIR)
    files = [f for f in files if f[-3:] == 'csv']
    files.sort()
    print(files)

    for f in files:
        print('Plotting %s' % f)
        full_date = f.split('_')[0]
        year = full_date[:4]
        if len(full_date) == 6:
            month = full_date[-2:]
            month_name = datetime.date(1900, int(month), 1).strftime('%B')
            title = '%s %s %s' (region_dict[region], month_name, year)
        else:
            title = '%s %s %s' (region_dict[region], full_date[:3:], year)

        fig, ax = fp.get_figax(maps=True, lats=latlim, lons=lonlim,
                               max_width=config.BASE_WIDTH,
                               max_height=config.BASE_HEIGHT)
        ax = fp.format_map(ax, latlim, lonlim, draw_labels=False)

        data = pd.read_csv(join(DATA_DIR, f))

        fig, ax, c = plot_TROPOMI(data, latlim, lonlim, res,
                                  figax=[fig, ax],
                                  vals='xch4',
                                    **plot_options)

        cax = fp.add_cax(fig, ax)
        cbar = fig.colorbar(c, cax=cax)
        cbar = fp.format_cbar(cbar, 'XCH4 (ppb)')
        fp.add_title(ax, title=title)

        # Save plot
        fp.save_fig(fig, DATA_DIR, '%s' % full_date)
