from os import listdir, mkdir, getcwd, environ
from os.path import join, dirname, realpath
import sys

import pandas as pd
import numpy as np
import xarray as xr
import math

import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from cartopy import feature as cf
sys.path.append('/n/home04/hnesser/TROPOMI/python')
import format_plots as fp
import config as config

# from fastkml import kml
# from lxml import etree
# from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
#                        AltitudeMode, Camera)

import cartopy.crs as ccrs
import cartopy

import shapefile as shp

# Tell matplotlib not to look for an X-window
environ['QT_QPA_PLATFORM']='offscreen'
rcParams['text.usetex'] = False

colors = plt.cm.get_cmap('inferno', lut=8)

# Other font details
rcParams['font.family'] = 'serif'
rcParams['font.size'] = config.LABEL_FONTSIZE*config.SCALE
#rcParams['text.usetex'] = True

def plot_TROPOMI(data, latlim, lonlim, res, figax=None, title='', genkml=False, vals='xch4', **plot_options):
    # Set default plot options]
    if 'vmin' not in plot_options:
        plot_options['vmin'] = 1700
    if 'vmax' not in plot_options:
        plot_options['vmax'] = 1900
    if 'cmap' not in plot_options:
        plot_options['cmap'] = 'inferno'

    res_prec = len(str(res).split('.')[1])

    # Define edges
    latlim = [round(l, res_prec) for l in latlim]
    lat_steps = int((latlim[1]-latlim[0])/res) + 1

    lonlim = [round(l, res_prec) for l in lonlim]
    lon_steps = int((lonlim[1]-lonlim[0])/res) + 1

    lats = np.around(np.linspace(latlim[0], latlim[1], lat_steps) + res/2,
                     res_prec + 1)
    lats_s = pd.DataFrame({'idx' : np.ones(len(lats)-1), 'lat' : lats[:-1]})
    lons = np.around(np.linspace(lonlim[0], lonlim[1], lon_steps) + res/2,
                     res_prec + 1)
    lons_s = pd.DataFrame({'idx' : np.ones(len(lons)-1), 'lon' : lons[:-1]})

    df = pd.merge(lats_s, lons_s, on='idx').drop(columns='idx')
    data = pd.merge(df, data, on=['lat', 'lon'], how='left')

    d_p = data.pivot(index='lat', columns='lon', values=vals)
    lon, lat = np.meshgrid(lons-res/2, lats-res/2)

    if figax is None:
        fig, ax = plt.subplots(figsize=(10,10),
                               subplot_kw={'projection' : ccrs.PlateCarree()})
    else:
        fig, ax = figax

    c = ax.pcolormesh(lons, lats, d_p,
                      vmin=plot_options['vmin'], vmax=plot_options['vmax'],
                      cmap=plot_options['cmap'], edgecolors=None, snap=True)

    return fig, ax, c

# def plot_TROPOMI(data, latlim, lonlim, res, figax=None, title='', genkml=False, vals='xch4', **plot_options):
#     # Set default plot options]
#     if 'vmin' not in plot_options:
#         plot_options['vmin'] = 1700
#     if 'vmax' not in plot_options:
#         plot_options['vmax'] = 1900
#     if 'cmap' not in plot_options:
#         plot_options['cmap'] = 'inferno'

#     # Define edges
#     lat_steps = (latlim[1]-latlim[0])/res + 1
#     lon_steps = (lonlim[1]-lonlim[0])/res + 1
#     if (lat_steps != int(lat_steps)) or (lon_steps != int(lon_steps)):
#         sys.exit('Bad interval.')

#     lats = np.around(np.linspace(latlim[0], latlim[1], int(lat_steps)) + res/2,
#                      len(str(res/2).split('.')[1]))
#     lats_s = pd.DataFrame({'idx' : np.ones(len(lats)-1), 'lat' : lats[:-1]})
#     lons = np.around(np.linspace(lonlim[0], lonlim[1], int(lon_steps)) + res/2,
#                      len(str(res/2).split('.')[1]))
#     lons_s = pd.DataFrame({'idx' : np.ones(len(lons)-1), 'lon' : lons[:-1]})
#     df = pd.merge(lats_s, lons_s, on='idx').drop(columns='idx')
#     data = pd.merge(df, data, on=['lat', 'lon'], how='left')

#     d_p = data.pivot(index='lat', columns='lon', values=vals)
#     lon, lat = np.meshgrid(lons-res/2, lats-res/2)

#     if figax is None:
#         fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection' : ccrs.PlateCarree()})
#     else:
#         fig, ax = figax
    
#     c = ax.pcolormesh(lons, lats, d_p,
#                       snap=True,
#                       vmin=plot_options['vmin'], vmax=plot_options['vmax'],
#                       cmap=plot_options['cmap'],
#                       edgecolors=None)

#     # Get kwargs
#     gridline_kwargs = {}
#     gridline_kwargs['draw_labels'] = gridline_kwargs.get('draw_labels', True)
#     gridline_kwargs['color'] = gridline_kwargs.get('color', 'grey')
#     gridline_kwargs['linestyle'] = gridline_kwargs.get('linestyle', ':')

#     # Format
#     gl = ax.gridlines(**gridline_kwargs)
#     gl.xlabel_style = {'fontsize' : config.TICK_FONTSIZE*config.SCALE}
#     gl.ylabel_style = {'fontsize' : config.TICK_FONTSIZE*config.SCALE}

#     return fig, ax, c

if __name__ == '__main__':
    DATA_DIR = sys.argv[1]
    REGION = str(sys.argv[2].split(',')[0])
    DATA_DIR = join(DATA_DIR, REGION)

    # DATA_DIR = '/n/seasasfs02/hnesser/TROPOMI/oversampling_output_csvs_14_14/dams/'
    # REGION = 'harsha'
    # DATA_DIR = join(DATA_DIR, REGION)

    region_dict = {'itaipu' : 'Itaipu',
                   'tucurui' : 'Tucurui',
                   'robertbourassa' : 'Robert-Bourassa (La Grande-2)',
                   'itumbiara' : 'Itumbiara',
                   'pampulha' : 'Pampulha',
                   'segredo' : 'Segredo',
                   'hartwell' : 'Hartwell',
                   'furnas' : 'Furnas',
                   'laforge' : 'Laforge-1',
                   'barrabonita' : 'Barra Bonita ',
                   'funil' : 'Funil',
                   'balbina' : 'Balbina',
                   'wattsbar' : 'Watts Bar',
                   'eastmain' : 'Eastmain 1',
                   'douglas' : 'Douglas',
                   'kariba' : 'Kariba',
                   'harsha' : 'Harsha',
                   'guntersville' : 'Guntersville',
                   'allatoona' : 'Allatoona',
                   'fontana' : 'Fontana'}

    #latlim = np.array(sys.argv[2].split(',')[1:3]).astype(float)
    #lonlim = np.array(sys.argv[2].split(',')[3:5]).astype(float)
    # lat = 39.02424
    # lon = -84.131436
    lat = float(sys.argv[2].split(',')[5])
    lon = float(sys.argv[2].split(',')[6])
    latlim = np.array([lat - 1, lat + 1]).round(2)
    lonlim = np.array([lon - 1, lon + 1]).round(2)
    res = 0.01
    count_min = 1

    files = listdir(DATA_DIR)
    files = [f for f in files if f[-3:] == 'csv']
    files.sort()

    for f in files:
        print(REGION)
        print(lat, lon)
        data = pd.read_csv(join(DATA_DIR, f))
        background = data['xch4'].quantile(0.1)
        data = data[(data['lon'] >= lonlim[0]) &
                    (data['lon'] <= lonlim[1]) &
                    (data['lat'] >= latlim[0]) &
                    (data['lat'] <= latlim[1])]
                
        # Calculate the number of data points expected
        max_count = ((latlim[1]-latlim[0])/res)*((lonlim[1]-lonlim[0])/res)
        real_count = data.shape[0]
        if real_count/max_count > 0.4:
            print('Plotting %s' % f)
            full_date = f.split('_')[0]
            year = full_date[:4]
            if len(full_date) == 6:
                month = full_date[-2:]
                month_name = datetime.date(1900, int(month), 1).strftime('%B')
                title = '%s %s %s' % (region_dict[REGION], month_name, year)
            else:
                title = '%s %s %s' % (region_dict[REGION], full_date[-3:], year)

            # Set up figure
            fig, ax = fp.get_figax(maps=True, lats=latlim, lons=lonlim,
                                   max_width=config.BASE_WIDTH,
                                   max_height=config.BASE_HEIGHT)
            ax = fp.format_map(ax, latlim, lonlim, draw_labels=True)
            cax = fp.add_cax(fig, ax)
            fp.add_title(ax, title=title, y=1.1)

            # Add features
            rivers = cf.NaturalEarthFeature('physical', 
                                            'rivers_lake_centerlines', 
                                            '10m')
            lakes = cf.NaturalEarthFeature('physical', 'lakes', '10m')
            ax.add_feature(rivers, facecolor='None', edgecolor='steelblue')
            ax.add_feature(lakes, facecolor='steelblue', edgecolor='steelblue')

            # Remove background
            data['enhancement'] = data['xch4'] - background

            # calculate vmin and vmax
            vmin = math.floor(data['xch4'].min()/50)*50
            vmax = math.ceil(data['xch4'].max()/50)*50
            #vmin = -np.max(np.abs([vmin, vmax]))
            #vmax = np.max(np.abs([vmin, vmax]))
            plot_options = {'vmin' : 1750,
                            'vmax' : 1950,
                            'cmap' : 'plasma'}

            fig, ax, c = plot_TROPOMI(data, latlim, lonlim, res,
                                      figax=[fig, ax], vals='xch4',
                                      **plot_options)
            ax.scatter(lon, lat, marker='x', s=50, color='black',
                       zorder=10)

            # Plot colorbar
            cbar = fig.colorbar(c, cax=cax)
            cbar = fp.format_cbar(cbar, 'XCH4 (ppb)')
            
            # Save plot
            fp.save_fig(fig, DATA_DIR, '%s_%s' % (REGION, full_date))
            plt.close()
            print('\n')
