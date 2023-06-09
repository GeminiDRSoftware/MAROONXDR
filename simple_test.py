import glob
import astrodata
import maroonx_instruments
from recipe_system.reduction.coreReduce import Reduce
from gempy.adlibrary import dataselect

# Get all files in the science_dir
all_files = glob.glob('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/science_dir/flats/*_r_0001.fits')

all_files.sort()
just_darks = dataselect.select_data(all_files, ['FLAT', 'DARK'])

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(all_files)
myreduce.drpkg= 'maroonxdr'
myreduce.runr()