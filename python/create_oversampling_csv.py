import xarray as xr
import magic
import pandas as pd
from os.path import join

def extract_tropomi_data(file, loc='./'):
    file_name = file.split('/')[-1][20:35]
    if magic.from_file(file) == 'Hierarchical Data Format (version 5) data':
        model_ds = xr.open_dataset(file, group='/PRODUCT')
        qa = model_ds[['qa_value']]
        if qa['qa_value'].where(qa['qa_value'] > 0.5, drop=True).shape == (0,0,0):
            print('No data with qa > 0.5.')

        else:
            model_ds = model_ds[['methane_mixing_ratio',\
                                 'methane_mixing_ratio_bias_corrected',\
                                 'methane_mixing_ratio_precision',\
                                 'qa_value']]

            g = 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS'
            model_ds = model_ds.merge(xr.open_dataset(file, group=g)[['surface_albedo_SWIR']])

            g = 'PRODUCT/SUPPORT_DATA/INPUT_DATA'
            model_ds = model_ds.merge(xr.open_dataset(file, group=g)[['surface_altitude']])
            model_ds = model_ds.where(model_ds['qa_value'] > 0.5, drop=True)

            for variable in model_ds.variables.values():
                if 'coordinates' in variable.attrs:
                    del variable.attrs['coordinates']

            model_ds.to_netcdf(loc + file_name + '.ch4')

            # corners
            model_ds = xr.open_dataset(file, group='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
            model_ds = model_ds.merge(qa)
            model_ds = model_ds.where(model_ds['qa_value'] > 0.5, drop=True)
            model_ds[['latitude_bounds',
                      'longitude_bounds']].to_netcdf(loc + file_name + '.corners')
    else:
        print('File is of unknown format.')

def make_csv(file, save_loc):
    file_name = file.split('/')[-1]
    if magic.from_file(file) == 'Hierarchical Data Format (version 5) data':
        # Start by getting the basic data: qa value, xch4, and xch4 uncertainty
        ds = xr.open_dataset(file, group='/PRODUCT')[['qa_value', 
                                                      'methane_mixing_ratio_bias_corrected',
                                                      'methane_mixing_ratio_precision']]
        
        # We want to avoid processing files where there are no data with qa > 0.5
        if ds['qa_value'].where(ds['qa_value'] > 0.5, drop=True).shape == (0,0,0):
            print('No data with qa > 0.5.')

        else:
            # We will merge in the other necessary data, namely the corners,
            # albedo, aod and number, etc.
            corners = xr.open_dataset(file,
                                      group='/PRODUCT/SUPPORT_DATA/GEOLOCATIONS')
            corners = corners[['latitude_bounds', 'longitude_bounds']]

            albedo = xr.open_dataset(file,
                                     group='PRODUCT/SUPPORT_DATA/DETAILED_RESULTS')
            albedo = albedo[['surface_albedo_SWIR',
                             'aerosol_optical_thickness_SWIR',
                             'aerosol_number_column']]

            pres = xr.open_dataset(file,
                                   group='/PRODUCT/SUPPORT_DATA/INPUT_DATA')
            pres = pres[['surface_pressure']]

            ds = ds.merge(corners)
            ds = ds.merge(albedo)
            ds = ds.merge(pres)

            # Get time information
            weekday = ds.time.dt.dayofweek.values[0]
            year = ds.time.dt.year.values[0]
            day = ds.time.dt.day.values[0]
            month = ds.time.dt.month.values[0]

            # Reformat the xarray 
            ds = ds.stack(z=['scanline', 'ground_pixel'])
            ds = ds.reset_index(dims_or_levels='z')
            ds = ds.squeeze('time').drop(['time', 'scanline', 'ground_pixel'])

            # Get rid of qa < 0.5
            ds = ds.where(ds['qa_value'] > 0.5, drop=True)

            # Get out corners and convert to dataframe
            latlons = ds[['latitude_bounds', 'longitude_bounds']].to_dataframe()

            lats = latlons.reset_index().pivot(index='z', columns='corner', values='latitude_bounds')
            lats.columns = ['lat' + str(col) for col in lats.columns]

            lons = latlons.reset_index().pivot(index='z', columns='corner', values='longitude_bounds')
            lons.columns = ['lon' + str(col) for col in lons.columns]

            # Now do the rest
            ds = ds.drop(['latitude_bounds', 'longitude_bounds'])
            ds = ds.to_dataframe()

            # Join in lat/lon corners
            ds = pd.concat([ds, lats, lons], axis=1)

            # Rename columns
            ds = ds.reset_index()
            ds = ds.rename(columns={'z' : 'idx',
                                    'qa_value' : 'qa',
                                    'methane_mixing_ratio_bias_corrected' : 'xch4',
                                    'methane_mixing_ratio_precision' : 'xch4_unc',
                                    'surface_albedo_SWIR' : 'albedo_swir',
                                    'aerosol_optical_thickness_SWIR' : 'aerosol_op',
                                    'aerosol_number_column' : 'aerosol_num',
                                    'surface_pressure' : 'surfpres',
                                    'latitude' : 'lat',
                                    'longitude' : 'lon'})

            # Insert columns for date, etc
            ds.insert(0, 'weekday', weekday)
            ds.insert(0, 'year', year)
            ds.insert(0, 'day', day)
            ds.insert(0, 'month', month)

            # Calculate ch4col
            ds['ch4col'] = (ds['xch4']*1e-9)*(ds['surfpres']/(28.97e-3*9.8))

            # Reorder the columns!
            ds = ds[['idx', 
                     'lat0', 'lat1', 'lat2', 'lat3', 'lat',
                     'lon0', 'lon1', 'lon2', 'lon3', 'lon',
                     'xch4', 'xch4_unc', 'qa',
                     'weekday', 'year', 'day', 'month',
                     'albedo_swir', 'aerosol_op', 'aerosol_num',
                     'ch4col', 'surfpres']]

            # Save a temporary csv
            temp_time = file_name.split('_')[8].split('T')[1]
            name = '%d%02d%02d_%s_temp.csv' % (year, month, day, temp_time) 
            ds.to_csv(join(save_loc, name), sep=',', header=True, index=False)

if __name__ == '__main__':
    from os import remove, listdir
    from os.path import join
    import sys

    data_dir = str(sys.argv[1])
    save_dir = str(sys.argv[2])
    min_date = int(sys.argv[3])

    files = listdir(data_dir)
    files = [f for f in files if (f[-2:] == 'nc')] 
    files = [f for f in files if (int(f[20:28]) >= min_date)]
    files.sort()

    ##start_num = int(sys.argv[1])
    ##end_num = int(sys.argv[2])

    ## Make the temporary csvs
    print('PROCESSING RAW TROPOMI INPUT.')
    for file in files:
        print('Processing ' + file.split('/')[-1])
        make_csv(join(data_dir, file), 
                 save_dir)

    # Now combine the csvs
    print('')
    print('COMBINING INPUT INTO MONTHLY CSVS')

    files = listdir(save_dir)
    files = [f for f in files if f[-8:] == 'temp.csv']
    files.sort()

    # Pull out the months of interest
    months = [f[:6] for f in files]
    months = list(set(months))
    months.sort()

    # Now combine the relevant months
    for m in months:
        print(m)
        mfiles = [f for f in files if f[:6] == m]
        mfiles.sort()
        mcsv = pd.concat([pd.read_csv(join(save_dir, f)) for f in mfiles])
        mcsv.to_csv(join(save_dir, '%s_combined.csv' % m), index=False, header=False)

        # Now create the shorter file that we will be using.
        mcsv = mcsv[['idx', 
                     'lat0', 'lat1', 'lat2', 'lat3', 'lat',
                     'lon0', 'lon1', 'lon2', 'lon3', 'lon',
                     'xch4', 'xch4_unc']]
        mcsv = mcsv[(mcsv['lat'] >= -80) & (mcsv['lat'] <= 80)]
        mcsv.to_csv(join(save_dir, '%s_combined_latlim.csv' % m), 
                    sep=',', 
                    float_format='%.6f',
                    header=False, 
                    index=False)
        for mf in mfiles:
            remove(join(save_dir, mf))
