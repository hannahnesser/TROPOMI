import xarray as xr
import magic

def extract_tropomi_data(file, loc='./'):
    file_name = file.split('/')[-1][0:35]
    if magic.from_file(file) == 'Hierarchical Data Format (version 5) data':
        model_ds = xr.open_dataset(file, group='/PRODUCT')
        qa = model_ds[['qa_value']]
        if qa['qa_value'].where(qa['qa_value'] > 0.5, drop=True).shape == (0,0,0):
            print('No data with qa > 0.5.')

        else:
            model_ds = model_ds[['methane_mixing_ratio',\
                                 'methane_mixing_ratio_bias_corrected',\
                                 'methane_mixing_ratio_precision',\
                                 'qa_value',\
                                 'delta_time']]

            g = 'PRODUCT/SUPPORT_DATA/DETAILED_RESULTS'
            model_ds = model_ds.merge(xr.open_dataset(file, group=g)[['surface_albedo_SWIR','column_averaging_kernel']])

            g = 'PRODUCT/SUPPORT_DATA/INPUT_DATA'
            model_ds = model_ds.merge(xr.open_dataset(file, group=g)[['surface_altitude','surface_pressure', 'pressure_interval','methane_profile_apriori']])
            model_ds = model_ds.where(model_ds['qa_value'] > 0.5, drop=True)

            for variable in model_ds.variables.values():
                if 'coordinates' in variable.attrs:
                    del variable.attrs['coordinates']



                    
            model_ds.to_netcdf(loc + file_name + '.nc')


    else:
        print('File is of unknown format.')

if __name__ == '__main__':
    from os import remove, listdir
    from os.path import join

    data_dir = '/n/seasasfs02/hnesser/TROPOMI/downloads_20190715/'
    save_dir = '/n/scratchlfs/jacob_lab/msulprizio/CH4/TROPOMI/processed_201908/'
    files = listdir(data_dir)
    files.sort()
    print(files[3000:])
    for file in files:
        print('Processing ' + file.split('/')[-1])
        extract_tropomi_data(join(data_dir, file), save_dir)
