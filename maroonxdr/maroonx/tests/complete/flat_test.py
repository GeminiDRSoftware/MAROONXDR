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

# from pathlib import Path

# from recipe_system.reduction.coreReduce import Reduce

# import maroonx_instruments  # noqa : important to load adclass tags

# # Get all files in the science_dir.  Change the path here to suit your installation.
# science_dir = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')

# FDDDF_files = list(science_dir.glob('*_DDDDF_r_0002.fits'))
# DFFFD_files = list(science_dir.glob('*_DFFFD_r_0002.fits'))

# FDDDF_files = [str(f) for f in FDDDF_files]
# DFFFD_files = [str(f) for f in DFFFD_files]

# FDDDF_files.sort()
# DFFFD_files.sort()

# # Run reduce on all selected files
# myreduce = Reduce()
# myreduce.files.extend(FDDDF_files)
# myreduce.files.extend(DFFFD_files)
# myreduce.drpkg = 'maroonxdr'
# myreduce.runr()


import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa : important to load adclass tags

logutils.config(file_name="test_reduction.log", stomp=False)
log = logutils.get_logger("test_reduction.log")
log.setLevel("DEBUG")


def test_reduce_dark():

    # Get all files in the science_dir.  Change the path here to suit your installation.
    test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
    science_dir = test_path / 'science_dir'

    # Get all files in the science_dir
    all_files = list(Path(science_dir).glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    # Change working directory to science_dir
    original_dir = os.getcwd()
    os.chdir(science_dir)

    for arm in ['BLUE', 'RED']:
        only_flats = dataselect.select_data(all_files, tags=['RAW', 'FLAT', arm])
        # FDDDF_files = list(science_dir.glob('*_DDDDF_r_0002.fits'))
        # DFFFD_files = list(science_dir.glob('*_DFFFD_r_0002.fits'))

        # FDDDF_files = [str(f) for f in FDDDF_files]
        # DFFFD_files = [str(f) for f in DFFFD_files]

        # FDDDF_files.sort()
        # DFFFD_files.sort()

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_flats)
        # myreduce.files.extend(FDDDF_files)
        # myreduce.files.extend(DFFFD_files)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeProcessedFlatDFFFF'
        myreduce.runr()

    # Restore original working directory
    os.chdir(original_dir)


if __name__ == '__main__':

    test_reduce_dark()
