from os import listdir, mkdir, getcwd, environ
from os.path import join, dirname, realpath
import sys

import pandas as pd
import numpy as np
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams

from fastkml import kml
from lxml import etree
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

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

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax

def plot_ong(data_dir, lat, lon, res=0.01,
             count_min=1, shape=None, **plot_options):
    files = listdir(data_dir)
    files = [f for f in files if (len(f.split('_')[0]) == 6)
                               & (f.split('.')[-1] == 'csv')]
    files.sort()
    
    for i, f in enumerate(files):
        fig, ax = plt.subplots(1, 2, figsize=((lat[1] - lat[0])*2*2 + 5, 
                                              (lon[1] - lon[0])*2),
                           subplot_kw={'projection' : ccrs.PlateCarree()})
        fig.suptitle('%s-%s' % (f.split('_')[0][:4], f.split('_')[0][-2:]),
                     fontsize=30)
        plt.subplots_adjust(wspace=0.1)
        
        data = pd.read_csv(join(data_dir, f + ''))
        data = data[data['cnt'] >= count_min]
        
        plot_options = {}
        fig, ax[0], c = plot_TROPOMI(data, lat, lon, res, 
                                     vals='xch4',
                                     figax=[fig, ax[0]],
                                     **plot_options)
        fig.colorbar(c, ax=ax[0], pad=0.1)
        
        _, deciles = pd.qcut(data['xch4'], q=10, retbins=True)
        background = deciles[1]
        data['enhancement'] = data['xch4'] - background
        plot_options['cmap'] = plot_options.pop('cmap', 'Reds')
        plot_options['vmin'] = 0
        plot_options['vmax'] = math.ceil(data['enhancement'].max()/10)*10
        fig, ax[1], c = plot_TROPOMI(data, lat, lon, res, 
                                     vals='enhancement',
                                     figax=[fig, ax[1]],
                                     **plot_options)
        fig.colorbar(c, ax=ax[1], pad=0.1)

        if shape is not None:
            for axis in ax:
                outline = plt.Polygon(p, fill=False, edgecolor='black', linewidth=2)
                axis.add_patch(outline)
        
        fig.savefig(join(data_dir, f.split('.')[0] + '.png'))

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

    # Get pngs
    files = listdir(data_dir)
    kml_pngs = [f for f in files if (f.split('.')[-1] == 'png')
                                  & (f != 'legend.png')
                                  & (f.split('.')[-2][-3:] != 'old')]
    files = [f for f in files if (len(f.split('_')[0]) == 7)
                               & (f.split('.')[-1] == 'csv')
                               & (f[:-4] + '.png' not in kml_pngs)]

    # Seasonal plotting
    #kml_pngs = []
    for i, f in enumerate(files):
        print('Plotting %s.' % f)
        data = pd.read_csv(join(data_dir, f))
        data = data[data['cnt'] >= count_min]
        
        plt.ioff()

        fig, ax = gearth_fig(llcrnrlon=latlon[2],
                             llcrnrlat=latlon[0],
                             urcrnrlon=latlon[3],
                             urcrnrlat=latlon[1],
                             pixels=1000)
        print('Basis figure generated.')

        fig, ax, c = plot_TROPOMI(data, 
                                  latlon[:2], 
                                  latlon[2:], 
                                  res,
                                  figax=[fig, ax],
                                  genkml=genkml,
                                  **plot_options)
        print('Figure complete.')

        ax.set_axis_off()
        
        save_name = f.split('.')[0] + '.png'
        kml_pngs.append(join(data_dir, save_name))
        fig.savefig(join(data_dir, save_name), 
                    transparent=True,
                    format='png')
        print('Figure saved at %s.' % join(data_dir, save_name))

        # Color bar time 
        if genkml and (i == (len(files) - 1)):
            fig = plt.figure(figsize=(1.5, 6.0), facecolor=None, frameon=False)
            ax = fig.add_axes([0.05, 0.05, 0.2, 0.9])
            cb = fig.colorbar(c, cax=ax)
            cb.set_label(r'XCH$_4$ [ppb]', rotation=-90, color='white', labelpad=20)
            cb.ax.tick_params(color='white', labelcolor='white')
            cb.outline.set_edgecolor('white')
            fig.savefig(join(data_dir, 'legend.png'), transparent=True, format='png')

    if genkml:
        # Now join them together as a kml
        print('Creating kmz from:')
        kml_pngs.sort()
        print(kml_pngs)
        make_kml(llcrnrlon=latlon[2], llcrnrlat=latlon[0],
                 urcrnrlon=latlon[3], urcrnrlat=latlon[1],
                 figs=kml_pngs, colorbar=join(data_dir, 'legend.png'),
                 kmzfile=join(data_dir, region + '_xch4_seasonal.kmz'), 
                 color='white')
