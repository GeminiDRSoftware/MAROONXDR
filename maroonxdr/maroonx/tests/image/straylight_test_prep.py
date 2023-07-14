import glob
import astrodata
from recipe_system.reduction.coreReduce import Reduce
from gempy.adlibrary import dataselect
import sys

from pathlib import Path
parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
from maroonxdr.maroonx.primitives_maroonx import MAROONX
import maroonx_instruments

"""
This is a file used to prepare straylight removed flats.  NOTE:  This does not produce a full flat.
Its purpose is to allow you to test the strayLight removal process.  Use this to create 2 flats and then take these
flats and place them in you science_dir folder.  Then run the test_stray_light_removal.py file to ensure everything is working as intended.
End users should not have to do this, unless something has gone wrong in the installation.
"""

# Get a DFFFD and FDDDF file from the science_dir.  Change the path here to suit your installation.
DFFFD_file = glob.glob('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/science_dir/20220725T162635Z_DFFFD_r_0001.fits')
FDDDF_file = glob.glob('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/science_dir/20220725T164012Z_FDDDF_r_0001.fits')

# Run reduce on all selected files to produce straylight_flats
myreduce = Reduce()
myreduce.files.extend(DFFFD_file)
myreduce.files.extend(FDDDF_file)
myreduce.drpkg= 'maroonxdr'
myreduce.recipename = 'makeStrayLightCheck'
myreduce.runr()