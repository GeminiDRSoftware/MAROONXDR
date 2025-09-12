import os
import itertools as it
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import astrodata
import maroonx_instruments  # noqa : important to load adclass tags


# =============================================================================
# Step 0 - Reduction setup
# =============================================================================
# Configure logging
logutils.config(file_name="test_reduction.log", mode="debug", stomp=True)

# Get all files in the science_dir. Change the path here to suit your installation.
science_dir = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')
os.chdir(science_dir)

# Calibration files paths
proc_dark = science_dir / "calibrations" / "processed_dark"
proc_flat = science_dir / "calibrations" / "processed_flat"

def get_files(path=None):
    """Get all files in a directory."""
    # If path is not provided, use the current working directory
    if path is None:
        path = science_dir
    all_files = list(path.glob('*.fits'))
    all_files = [str(f) for f in all_files]
    all_files.sort()
    return all_files

# Define tags for file selection
exptime_tags = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
arm_tags = ['BLUE', 'RED']


# =============================================================================
# Step 1 - Debundle the data
# =============================================================================
# Select bundles
selected_bundles = dataselect.select_data(get_files(), tags=['BUNDLE'])

# Run reduce on all selected files
# Dragons will choose the correct recipe based on the tags
myreduce = Reduce()
myreduce.files.extend(selected_bundles)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()


# =============================================================================
# Step 2 - Master Flats
# =============================================================================
# Flats should be run for red and blue arms separately
for arm in arm_tags:

    # Select both DFFFD and FDDDF files
    selected_flats = dataselect.select_data(get_files(), tags=['RAW', 'FLAT', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_flats)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeProcessedFlatDFFFF'
    myreduce.runr()


# =============================================================================
# Step 3 - Master Darks
# =============================================================================
# Dark frames should be run for each exposure time and arm
for exptime, arm in it.product(exptime_tags, arm_tags):

    selected_darks = dataselect.select_data(get_files(), tags=['RAW', 'DARK', exptime, arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()


# =============================================================================
# Step 4 - Create Coefficient Darks
# =============================================================================
for arm in arm_tags:

    # Select master darks with PROCESSED tag 
    selected_darks = dataselect.select_data(
        get_files(proc_dark), tags=['PROCESSED', 'DARK', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeDarkCoefficients'
    myreduce.runr()


# =============================================================================
# Step 5 - Synthetic Darks for science frames
# =============================================================================
# synthetic darks
for arm in arm_tags:

    selected_darks = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_darks)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeSyntheticDark'
    myreduce.runr()


# =============================================================================
# Step 6 - Wavecal frames
# =============================================================================
for arm in arm_tags:

    selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'WAVECAL', arm])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_spect)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'makeDynamicWavecal'
    myreduce.runr()


# =============================================================================
# Step 7 - Science frames
# =============================================================================
for arm in arm_tags:

    selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm, '300s'])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_spect)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()
