from pathlib import Path
import itertools as it

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import astrodata
import maroonx_instruments  # noqa : important to load adclass tags

# Configure logging
logutils.config(file_name="test_reduction.log", mode="debug", stomp=True)

# Get all files in the science_dir.  Change the path here to suit your installation.
science_dir = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')

def get_files(path=None):
    if path is None:
        path = science_dir
    all_files = list(path.glob('*.fits'))
    all_files = [str(f) for f in all_files]
    all_files.sort()
    return all_files

# some testing files from previous runs are in science_dir
# so we need to filter them out, use RAW tag to select only raw files
raw_files = dataselect.select_data(get_files(), tags=['RAW'])


# =============================================================================
# Step 1 - Debundle the data
# Select bundles
selected_bundles = dataselect.select_data(get_files(), tags=['BUNDLE'])

# Run reduce on all selected files
# Dragons will choose the correct recipe based on the tags
myreduce = Reduce()
myreduce.files.extend(selected_bundles)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()

# Step 2 - Master Flats
# Flats should be run for red and blue arms separately
for arm in ['RED', 'BLUE']:

    # Select both DFFFD and FDDDF files
    selected_flats = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_flats)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeProcessedFlatDFFFF'
    myreduce.runr()


# Step 3 - Master Darks
# Dark frames should be run for each exposure time and arm
# exptime_tags = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
exptime_tags = ['300s']
arm_tags = ['BLUE', 'RED']

for exptime, arm in it.product(exptime_tags, arm_tags):

    selected_darks = dataselect.select_data(get_files(), tags=['RAW', 'DARK', exptime, arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()


# Step 4 - Extract flux
arm_tags = ['RED', 'BLUE']
for arm in arm_tags:

    selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'WAVECAL', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_spect)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeDynamicWavecal'
    myreduce.runr()
