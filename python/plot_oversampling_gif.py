import imageio
from os import listdir
from os.path import join

##########################
# Set default file paths
BASE_DIR = '/Users/hannahnesser/Documents/Harvard/Research/TROPOMI/oversampling_output/new_data'
DATA_DIR = join(BASE_DIR, 'northamerica')
PLOT_DIR = join(BASE_DIR, 'plots')
##########################

files = ['2018%02d' % i for i in range(5, 13)]
files = files + ['2019%02d' % i for i in range(1, 5)]
print(files)
files = ['TROPOMI_%s.png' % f for f in files]
files.sort()
print(files)

images = []
for f in files:
    images.append(imageio.imread(join(PLOT_DIR, f)))
imageio.mimsave(join(PLOT_DIR, 'northamerica.gif'), images,
                **{'duration' : 1})
