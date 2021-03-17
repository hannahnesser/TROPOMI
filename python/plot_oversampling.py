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

# Tell matplotlib not to look for an X-window
environ['QT_QPA_PLATFORM']='offscreen'

colors = plt.cm.get_cmap('inferno', lut=8)

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 14

def plot_TROPOMI(data, latlim, lonlim, res, figax=None, title='', genkml=False, vals='xch4', **plot_options):
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

    # Set default plot options
    plot_options['vmin'] = plot_options.pop('vmin', math.floor(data[vals].min()/50)*50)
    plot_options['vmax'] = plot_options.pop('vmax', math.ceil(data[vals].max()/50)*50)
    plot_options['cmap'] = plot_options.pop('cmap', 'inferno')

    if figax is None:
        fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection' : ccrs.PlateCarree()})
    else:
        fig, ax = figax

    c = ax.pcolormesh(lons, lats, d_p,
                      snap=True,
                      vmin=plot_options['vmin'], vmax=plot_options['vmax'],
                      cmap=plot_options['cmap'],
                      edgecolors=None)

    if not genkml:
        ax.set_title(title, fontsize=15, y=1.05)
        ax.set_extent(lonlim + latlim)
        ax.add_wms(wms='http://vmap0.tiles.osgeo.org/wms/vmap0',
                   layers=['basic'])
        ax.gridlines(linestyle=':', draw_labels=True, color='grey')

    return fig, ax, c

def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
#         draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 0)
        ground.name = kw.pop('name', fig.split('/')[-1].split('.')[0])
        ground.color = kw.pop('color', 'white')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.west = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.east = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='Color Bar')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=1, y=1,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.95, y=0.95,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)

# def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
#     """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
#     aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
#     xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
#     ysize = np.ptp([urcrnrlat, llcrnrlat])
#     aspect = ysize / xsize

#     if aspect > 1.0:
#         figsize = (10.0 / aspect, 10.0)
#     else:
#         figsize = (10.0, 10.0 * aspect)

#     if False:
#         plt.ioff()  # Make `True` to prevent the KML components from poping-up.
#     fig = plt.figure(figsize=figsize,
#                      frameon=False,
#                      dpi=pixels//10)
#     # KML friendly image.  If using basemap try: `fix_aspect=False`.
#     ax = fig.add_axes([0, 0, 1, 1])
#     ax.set_xlim(llcrnrlon, urcrnrlon)
#     ax.set_ylim(llcrnrlat, urcrnrlat)
#     return fig, ax

if __name__ == '__main__':
    data_dir = sys.argv[1]
    region   = str(sys.argv[2].split(',')[0])
    data_dir = join(data_dir, region)
    print('Processing %s' % region)

    if data_dir.split('/')[-2] == 'world':
        genkml = False
        print('Generating Google Earth Files.')
    else:
        genkml = False

    latlon = np.array(sys.argv[2].split(',')[1:]).astype(float)

    res  = 0.01
    vmin = 1700
    vmax = 1900
    plot_options = {'vmin' : vmin, 'vmax' : vmax}
    count_min = 1

    files = listdir(DATA_DIR)
    files = [f for f in files if f[-3:] == 'csv']
    files.sort()

    print(files)

    # for f in files:
    #     data = pd.read_csv(join(DATA_DIR, f))
    #     data = data[(data['lon'] >= lonlim[0]) &
    #                 (data['lon'] <= lonlim[1]) &
    #                 (data['lat'] >= latlim[0]) &
    #                 (data['lat'] <= latlim[1])]

    #     print('Plotting %s' % f)
    #     full_date = f.split('_')[0]
    #     year = full_date[:4]
    #     if len(full_date) == 6:
    #         month = full_date[-2:]
    #         month_name = datetime.date(1900, int(month), 1).strftime('%B')
    #         title = '%s %s %s' % (region_dict[REGION], month_name, year)
    #     else:
    #         title = '%s %s %s' % (region_dict[REGION], full_date[-3:], year)

    #     fig, ax = fp.get_figax(maps=True, lats=latlim, lons=lonlim,
    #                            max_width=config.BASE_WIDTH,
    #                            max_height=config.BASE_HEIGHT)
    #     ax = fp.format_map(ax, latlim, lonlim, draw_labels=False)

    #     # calculate vmin and vmax
    #     vmin = 1800
    #     vmax = 1900
    #     plot_options = {'vmin' : vmin,
    #                     'vmax' : vmax,
    #                     'cmap' : 'plasma'}

    #     fig, ax, c = plot_TROPOMI(data, latlim, lonlim, res,
    #                               figax=[fig, ax],
    #                               vals='xch4',
    #                               **plot_options)

    #     cax = fp.add_cax(fig, ax)
    #     cbar = fig.colorbar(c, cax=cax)
    #     cbar = fp.format_cbar(cbar, 'XCH4 (ppb)')
    #     fp.add_title(ax, title=title, y=1.1)

    #     # Save plot
    #     fp.save_fig(fig, DATA_DIR, '%s_%s' % (REGION, full_date))
    #     plt.close()
    #     print('\n')
