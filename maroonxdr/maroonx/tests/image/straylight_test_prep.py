"""Prepare straylight removed flats for testing.

This script prepares straylight removed flats for testing the straylight removal
process. It does not produce a full flat. Use this to create 2 flats and place them
in your science_dir folder, then run test_stray_light_removal.py to verify the process.
"""

from pathlib import Path

from recipe_system.reduction.coreReduce import Reduce

# Test data should be under science_dir
science_dir = Path(__file__).parents[4] / 'science_dir'

# Get all files in the science_dir
FDDDF_file = [str(science_dir / "20241114T190714Z_DDDDF_b_0007.fits")]
DFFFD_file = [str(science_dir / "20241114T182328Z_DFFFD_b_0008.fits")]

# Run reduce on all selected files to produce straylight_flats
myreduce = Reduce()
myreduce.files.extend(DFFFD_file)
myreduce.files.extend(FDDDF_file)
myreduce.drpkg= 'maroonxdr'
myreduce.recipename = 'makeStrayLightCheck'
myreduce.runr()
