"""
This script tests the creation of 1D reduced spectra for MAROON-X DEEEE fits files.  
It does not rely on pytest, and does not produce a success or fail output like pytest does.
Instead, if the reduce runs successfully, this will produce a 1D reduced spectra file that can be 
used for dynamic wavelength calibration.  End users should use this test to test their installation.
"""

import glob
import sys
import astrodata
from recipe_system.reduction.coreReduce import Reduce
from gempy.adlibrary import dataselect

from pathlib import Path
parent_dir = Path(__file__).parents[4]
sys.path.append(str(parent_dir))
import maroonx_instruments

all_files = glob.glob('/Users/rohan/Desktop/mxcrunch/outputs/*_DEEEE_r_*.fits')
    
all_files.sort()

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(all_files)
myreduce.drpkg= 'maroonxdr'
myreduce.runr()