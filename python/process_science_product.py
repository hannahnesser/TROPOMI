import xarray as xr
import pandas as pd
from datetime import date
from os import listdir
from os.path import join

# Create a function to process each file 
def process_science_product(file_name,
                            raw_data_dir='/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/raw',
                            processed_data_dir='/n/holyscratch01/jacob_lab/hnesser/TROPOMI/downloads_new/processed',
                            today=date.today().strftime('%Y%m%dT%H%M%S')):
    # Iterate through those files.
    print('Processing %s' % file_name)
    file = join(raw_data_dir, file_name)
    data = {}

    # First, open the diagnostics group and retrieve the qa_value
    # and a mask where the qa_value is not nan. We will apply
    # this mask to all future variables.
    d = xr.open_dataset(file, group='diagnostics')
    mask_qa = (d['qa_value'] <= 1)

    if mask_qa.sum() > 0:
        data['qa_value'] = d['qa_value'].where(mask_qa, drop=True)
        d.close()

        # Second, open the instrument group and retrieve information
        # on the time, lat/lon coordinates, and glint status.
        d = xr.open_dataset(file, group='instrument')
        vars = ['time', 'latitude_center', 'longitude_center', 'latitude_corners',
                'longitude_corners', 'glintflag']
        for var in vars:
            data[var] = d[var].where(mask_qa, drop=True)
            d.close()

        # Third, open the meteo group and retrieve information on the
        # altitude/pressure.
        d = xr.open_dataset(file, group='meteo')
        vars = ['altitude_levels', 'surface_altitude', 'dp', 'surface_pressure',
                'dry_air_subcolumns']
        for var in vars:
            data[var] = d[var].where(mask_qa, drop=True)
            d.close()

        # Fourth, open the target product group and retrieve all variables.
        # This includes xch4, xch4_precision, xch4_column_averaging_kernel,
        # ch4_profile_apriori, xch4_apriori, and xch4_corrected.
        d = xr.open_dataset(file, group='target_product')
        for var in d.keys():
            data[var] = d[var].where(mask_qa, drop=True)
            d.close()

        # Fifth, open the side product group and retrieve albedo and AOD
        d = xr.open_dataset(file, group='side_product')
        vars = ['surface_albedo', 'aerosol_optical_thickness']
        for var in vars:
            data[var] = d[var].where(mask_qa, drop=True)
            d.close()

        # Now that we've gotten all the data, we need to get information
        # on the date for the file name.
        dates = pd.DataFrame(data['time'].values[:,:-1],
                             columns=['year', 'month', 'day', 'hour', 'minute', 'second'])
        dates = pd.to_datetime(dates).dt.strftime('%Y%m%dT%H%M%S')
        start_date = dates.min()
        end_date = dates.max()

        # And we need the orbit number
        orbit_number = file_name.split('_')[-1].split('.')[0]

        # Now save the file name, using 0s instead of the collection number
        # and processor number. We also use ACMG as the processing stream (ha)
        # and the time of processing as the current time.
        save_name = ('S5P_ACMG_L2__CH4____%s_%s_%s_00_000000_%s.nc'
                     % (start_date, end_date, orbit_number, today))

        # Finally, merge data and save out.
        data = xr.merge(data.values())
        data.to_netcdf(join(processed_data_dir, save_name))

if __name__ == "__main__":
    import sys
    raw_data_dir = sys.argv[1]
    processed_data_dir = sys.argv[2]
    today = sys.argv[3]
    files = sys.argv[4:]
    for f in files:
        process_science_product(f,
                                raw_data_dir=raw_data_dir,
                                processed_data_dir=processed_data_dir,
                                today=today)

                                
