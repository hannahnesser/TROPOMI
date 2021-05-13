import xarray as xr
import pandas as pd
from datetime import date
import re
from os import listdir
from os.path import join

# Create a function to process each file
def process_science_product(file_name,
                            raw_data_dir='/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/raw',
                            processed_data_dir='/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/processed',
                            today=date.today().strftime('%Y%m%dT%H%M%S')):
    # Create a list of existing orbits
    preprocessed = [re.split('_+', f)[6] for f in listdir(processed_data_dir)]

    # Iterate through those files.
    # print('Processing %s' % file_name)
    file = join(raw_data_dir, file_name)
    data = {}

    # We need the orbit number
    orbit_number = file_name.split('_')[-1].split('.')[0]

    # First, open the diagnostics group and retrieve the qa_value
    # and a mask where the qa_value is not nan. We will apply
    # this mask to all future variables.
    d = xr.open_dataset(file, group='diagnostics')
    mask_qa = (d['qa_value'] <= 1)

    # Now get a mask for the masked values
    d = xr.open_dataset(file, group='target_product')
    mask_mk = (d['xch4_corrected'] != 9.96921e36)

    # Combine the masks
    mask = (mask_qa and mask_mk)

    # Where the qa value is <= 1 is a good measure of the retrieval success.
    total_nobs = len(mask)
    success_nobs = mask.sum().values

    if (mask.sum() > 0) and (orbit_number not in preprocessed):
        data['qa_value'] = d['qa_value'].where(mask, drop=True)
        d.close()

        # Second, open the instrument group and retrieve information
        # on the time, lat/lon coordinates, and glint status.
        d = xr.open_dataset(file, group='instrument')
        vars = ['time', 'latitude_center', 'longitude_center',
                'latitude_corners', 'longitude_corners', 'glintflag']
        for var in vars:
            data[var] = d[var].where(mask, drop=True)
            d.close()

        # Third, open the meteo group and retrieve information on the
        # altitude/pressure as well as cloud fraction
        d = xr.open_dataset(file, group='meteo')
        vars = ['altitude_levels', 'surface_altitude', 'dp',
                'surface_pressure', 'dry_air_subcolumns', 'cloud_fraction']
        for var in vars:
            data[var] = d[var].where(mask, drop=True)
            d.close()

        # Fourth, open the target product group and retrieve all variables.
        # This includes xch4, xch4_precision, xch4_column_averaging_kernel,
        # ch4_profile_apriori, xch4_apriori, and xch4_corrected.
        d = xr.open_dataset(file, group='target_product')
        for var in d.keys():
            data[var] = d[var].where(mask, drop=True)
            d.close()

        # Fifth, open the side product group and retrieve albedo and AOD
        d = xr.open_dataset(file, group='side_product')
        vars = ['surface_albedo', 'aerosol_optical_thickness']
        for var in vars:
            data[var] = d[var].where(mask, drop=True)
            d.close()

        # Now that we've gotten all the data, we need to get information
        # on the date for the file name.
        dates = pd.DataFrame(data['time'].values[:,:-1],
                             columns=['year', 'month', 'day', 'hour',
                                      'minute', 'second'])
        dates = pd.to_datetime(dates).dt.strftime('%Y%m%dT%H%M%S')
        start_date = dates.min()
        end_date = dates.max()

        # Now save the file name, using 0s instead of the collection number
        # and processor number. We also use ACMG as the processing stream (ha)
        # and the time of processing as the current time.
        save_name = ('S5P_ACMG_L2__CH4____%s_%s_%s_00_000000_%s.nc'
                     % (start_date, end_date, orbit_number, today))

        # Finally, merge data and save out.
        data = xr.merge(data.values())
        data.to_netcdf(join(processed_data_dir, save_name))

    return orbit_number, total_nobs, success_nobs

if __name__ == "__main__":
    import sys
    raw_data_dir = sys.argv[1]
    proc_data_dir = sys.argv[2]
    today = sys.argv[3]
    files = sys.argv[4:]
    print('ORBIT, NOBS_TOT, NOBS_SUCCESS   ')
    for f in files:
        o, t, s = process_science_product(f,
                                          raw_data_dir=raw_data_dir,
                                          processed_data_dir=proc_data_dir,
                                          today=today)
        print(f'{o}, {t}, {s}')

