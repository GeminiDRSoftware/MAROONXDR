"""
Script used to test the creation of flats for MAROON-X data.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a calibrated dark
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.

To run this test, simply run the command "python flat_test.py" in the terminal, after
ensuring the correct path to science_dir in Path().  This test expects you to have a
science_dir in the root directory of the installation (you have to make this directory)
and to have got the test fits files from Kathleen.
"""

from pathlib import Path

from recipe_system.reduction.coreReduce import Reduce

# Get all files in the science_dir.  Change the path here to suit your installation.
science_dir = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')

FDDDF_files = list(science_dir.glob('*_DDDDF_r_0002.fits'))
DFFFD_files = list(science_dir.glob('*_DFFFD_r_0002.fits'))

FDDDF_files.sort()
DFFFD_files.sort()

# Run reduce on all selected files
myreduce = Reduce()
myreduce.files.extend(FDDDF_files)
myreduce.files.extend(DFFFD_files)
myreduce.drpkg = 'maroonxdr'
myreduce.runr()
