import xarray as xr
import magic
import pandas as pd
import re
from os.path import join

def process_TROPOMI(data):
    # Filter on qa_value and masked methane values
    data = data.where(data['qa_value'] > 0.5, drop=True)
    data = data.where(data['xch4_corrected'] != 9.96921e36, drop=True)

    # Select only relevant variables
    # data = data[['time', 'qa_value',
    #              'latitude_center', 'longitude_center',
    #              'latitude_corners', 'longitude_corners',
    #              'xch4_corrected', 'xch4_precision',
    #              'surface_albedo', 'aerosol_optical_thickness',
    #              'surface_pressure']]
    data = data[['latitude_center', 'longitude_center',
                 'latitude_corners', 'longitude_corners',
                 'xch4_corrected', 'xch4_precision']]

    # Albedo and AOD have two columns [NIR, SWIR]. We select SWIR.
    # data = data.where(data.nwin == 1, drop=True).squeeze()

    # Get date information (this function used to get day of week,
    # we are discarding it)
    # data = data.assign(year=data['time'][:, 0],
    #                    month=data['time'][:, 1],
    #                    day=data['time'][:, 2])
    # data = data.drop('time')

    # And fix the latitude/longitude corner situation
    data = data.assign(lat0=data['latitude_corners'][:, 0],
                       lat1=data['latitude_corners'][:, 1],
                       lat2=data['latitude_corners'][:, 2],
                       lat3=data['latitude_corners'][:, 3],
                       lon0=data['longitude_corners'][:, 0],
                       lon1=data['longitude_corners'][:, 1],
                       lon2=data['longitude_corners'][:, 2],
                       lon3=data['longitude_corners'][:, 3])
    data = data.drop(labels=['latitude_corners', 'longitude_corners'])

    # # Calculate the column
    # data = data.assign(col=(data['xch4_corrected']
    #                         *data['surface_pressure']
    #                         *1e-9/28.97e-3*9.8))

    # Rename variables
    # data = data.rename({'nobs' : 'idx',
                        # 'qa_value' : 'qa',
    data = data.rename({'latitude_center' : 'lat',
                        'longitude_center' : 'lon',
                        'xch4_corrected' : 'xch4',
                        'xch4_precision' : 'xch4_unc'})
                        # 'surface_albedo' : 'albedo',
                        # 'aerosol_optical_thickness' : 'aod',
                        # 'surface_pressure' : 'pres'})

    return data

def make_csv(file, save_loc, save_name):
    file_name = file.split('/')[-1]
    data = data.to_dataframe().reset_index()

    # Reorder the columns!
    data = data[['nobs',
                 'lat0', 'lat1', 'lat2', 'lat3', 'lat',
                 'lon0', 'lon1', 'lon2', 'lon3', 'lon',
                 'xch4', 'xch4_unc', 'qa',
                 'year', 'day', 'month',
                 'albedo', 'aod', 'col', 'pres']]

    # Save a temporary csv
    ds.to_csv(join(save_loc, save_name), sep=',', header=True, index=False)

if __name__ == '__main__':
    from os import remove, listdir
    from os.path import join
    import sys

    data_dir = str(sys.argv[1])
    save_dir = str(sys.argv[2])
    min_date = int(sys.argv[3])

    files = listdir(data_dir)
    files = [f for f in files if (f[-2:] == 'nc')]
    files = [join(data_dir, f) for f in files if (int(f[20:28]) >= min_date)]
    files.sort()

    # Separate into monthly
    file_dict = {}
    for file_name in files:
        # Get the date (YYYY, MM, and DD) of the raw TROPOMI file
        shortname = re.split('\/|\.', file_name)[-2]
        strdate = re.split('_+|T', shortname)
        start_month = strdate[4][:6]
        end_month = strdate[6][:6]

        # Add the file to the list of Sat_files
        if start_month in file_dict.keys():
            file_dict[start_month].append(file_name)
        else:
            file_dict[start_month] = [file_name]

        if start_month != end_month:
            if end_month in file_dict.keys():
                file_dict[end_month].append(file_name)
            else:
                file_dict[end_month] = [file_name]

    ##start_num = int(sys.argv[1])
    ##end_num = int(sys.argv[2])
    for month, files in file_dict.items():
        print('Processing %s' % month)
        ds = xr.open_mfdataset(files, concat_dim='nobs',
                               preprocess=process_TROPOMI,
                               combine='nested')
        ds = ds.to_dataframe().reset_index()

        # Reorder the columns!
        # ds = ds[['idx',
        #          'lat0', 'lat1', 'lat2', 'lat3', 'lat',
        #          'lon0', 'lon1', 'lon2', 'lon3', 'lon',
        #          'xch4', 'xch4_unc', 'qa',
        #          'year', 'day', 'month',
        #          'albedo', 'aod', 'col', 'pres']]

        # Save csv
        # ds.to_csv(join(save_loc, '%s.csv' % month),
        #           sep=',', header=True, index=False)
        ds = ds[['idx',
                 'lat0', 'lat1', 'lat2', 'lat3', 'lat',
                 'lon0', 'lon1', 'lon2', 'lon3', 'lon',
                 'xch4', 'xch4_unc']]
        ds = ds[(ds['lat'] >= -80) & (ds['lat'] <= 80)]
        ds.to_csv(join(save_dir, '%s_latlim.csv' % m),
                  sep=',',
                  float_format='%.6f',
                  header=False,
                  index=False)
