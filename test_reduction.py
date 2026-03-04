import os
import itertools as it
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import astrodata
import maroonx_instruments  # noqa : important to load adclass tags

# Run with Legacy patch?
LEGACY = False

# Get all files in the science_dir. Change the path here to suit your installation.
if LEGACY:
    science_dir = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir_legacy_patch')
else:
    science_dir = Path('/home/martin/Projects/MaroonX/MAROONXDR/science_dir')
os.chdir(science_dir)

# =============================================================================
# Step 0 - Reduction setup
# =============================================================================
# Configure logging
if LEGACY:
    logutils.config(file_name="test_reduction_legacy.log", mode="debug", stomp=True)
else:
    logutils.config(file_name="test_reduction.log", mode="debug", stomp=True)

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
    myreduce.uparms = {
        'removeStrayLight:legacy': LEGACY,
    }
    myreduce.runr()

    # Select masterflats
    selected_mflats = dataselect.select_data(get_files(), tags=['PROCESSED', 'FLAT', arm])

    myreduce = Reduce()
    myreduce.files.extend(selected_mflats)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'measureBlaze'
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
    # Exclude already created darks with DARK_SYNTH or DARK_COEFF tags
    selected_darks = dataselect.select_data(
        get_files(proc_dark), 
        tags=['PROCESSED', 'DARK', arm],
        xtags=['DARK_COEFF', 'DARK_SYNTH'])

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
sci_exptime_tags = ['120s', '300s', '600s', '900s', '1800s']

# for arm in arm_tags:
#     selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm, '300s'])
for exptime, arm in it.product(sci_exptime_tags, arm_tags):
    selected_spect = dataselect.select_data(get_files(), tags=['RAW', 'SCI', arm, exptime])

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(selected_spect)
    myreduce.drpkg = 'maroonxdr'
    myreduce.uparms = {
        'extractStripes:legacy': LEGACY,
    }
    myreduce.runr()

# =============================================================================
# Step 8 - Export Science Bundles
# =============================================================================

selected_spect = dataselect.select_data(
    get_files(), 
    tags=['PROCESSED', 'SCI', '300s'])

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(selected_spect)
myreduce.drpkg = 'maroonxdr'
myreduce.recipename = 'exportReducedBundle'
myreduce.uparms = {
    'bundleArmStreams:suffix': '_exported',
}
myreduce.runr()


# =============================================================================
# Step 9 - Barycentric Velocity Correction (with target-specific parameters)
# =============================================================================

# Define target-specific barycentric correction parameters
BARYCOR_TARGETS = [
    # (target_name, simbad_target_name, use_coords, zp_frd, zp_pc, exptime_tag)
    ('HD 203030', 'HD 203030', False, 0., 0., '300s'),
    ('Ross 248', 'Ross 248', False, 0., 0., '1800s'),
    ('Gl 12', 'G 32-5', False, 0., 0., '1800s'),
    ('TOI 4527', 'G 2-33', False, 0., 0., '1800s'),
    ('TOI 3686', 'TYC 3699-237-1', False, 0., 0., '600s'),
    ('TOI 4529', None, True, 0., 0., '1800s'),  # use_coords when SIMBAD fails
    ('NGTS-11', 'NGTS-11', False, 0., 0., '900s'),
    ('TOI-6662.01', None, True, 1.0, 7.7, '900s'),  # custom zeropoints
    ('TOI 4353', 'CD-34 1169', False, 0., 0., '1800s'),
    ('TOI-1634', 'PM J03455+3706', False, 0., 0., '600s'),
    ('TOI 4324', None, True, 0., 0., '1800s'),  # use_coords
    ('Gaia DR2 1146658793750549248', 'LSPM J1006+8305', False, 0., 0., '1800s'),
    ('Lalande 21185', 'HD 95735', False, 0., 0., '120s'),
]

for target_name, simbad_name, use_coords, zp_frd, zp_pc, exptime in BARYCOR_TARGETS:
    # Select reduced science files for this target and exposure time
    selected_files = dataselect.select_data(
        get_files(),
        tags=['PROCESSED', 'SCI', exptime]
    )
    # Build uparms dict for this target
    uparms = {
        'barycentricCorrection:target_name': target_name,
        'barycentricCorrection:use_coords': use_coords,
        'barycentricCorrection:simbad_target_name': simbad_name,
        'barycentricCorrection:zp_frd': zp_frd,
        'barycentricCorrection:zp_pc': zp_pc,
        'barycentricCorrection:report': True,
    }
    myreduce = Reduce()
    myreduce.files.extend(selected_files)
    myreduce.drpkg = 'maroonxdr'
    myreduce.recipename = 'applyBarycentricCorrection'
    myreduce.uparms = uparms
    myreduce.runr()