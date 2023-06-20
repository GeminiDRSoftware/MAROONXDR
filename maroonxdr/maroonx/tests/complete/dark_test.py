"""
This script is used to test the creation of darks for MAROON-X data.  It does not rely on pytest, and does not
produce a success or fail output like pytest does.  Instead, if the reduce runs successfully, this will produce
a calibrated dark file.  End users should use this test to test their installation. Only if there is an error
in this test (or other "complete" tests) should they use the echelle and image unit tests to test their installation.
To run this test, simply run the command "python dark_test.py" in the terminal, after ensuring the correct path in line
22.  This test expects you to have a science_dir in the root directory of the installation (you have to make this directory)
and to have got the test fits files from Kathleen.
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

# Get all files in the science_dir.  Change the path here to suit your installation.
all_files = glob.glob('/Users/rohan/Desktop/DRAGONS-X/MAROONXDR/science_dir/*_r_0300.fits')

all_files.sort()
just_darks = dataselect.select_data(all_files, ['DARK'])

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(just_darks)
myreduce.drpkg= 'maroonxdr'
myreduce.runr()