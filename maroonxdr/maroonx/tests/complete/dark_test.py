"""
Script used to test the creation of darks for MAROON-X data.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a calibrated dark
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.

To run this test, simply run the command "python dark_test.py" in the terminal, after
ensuring the correct path to science_dir in Path().  This test expects you to have a
science_dir in the root directory of the installation (you have to make this directory)
and to have got the test fits files from Kathleen.
"""

# from pathlib import Path

# from gempy.adlibrary import dataselect
# from recipe_system.reduction.coreReduce import Reduce

# import maroonx_instruments  # noqa : important to load adclass tags

# # Get all files in the science_dir.  Change the path here to suit your installation.
# science_dir = Path('/home/martin/Documentos/Projects/MAROONXDR/science_dir')
# all_files = list(science_dir.glob('*_DDDDE_*.fits'))
# all_files = [str(f) for f in all_files]
# all_files.sort()

# for arm_tag in ['BLUE', 'RED']:
#     just_darks = dataselect.select_data(all_files, tags=['DARK', arm_tag, '300s'])
#     print('Number of DARK files:', len(just_darks))

#     # Run reduce on all selected files
#     myreduce = Reduce()
#     myreduce.files.extend(just_darks)
#     myreduce.drpkg = 'maroonxdr'
#     # coment out this line for default reduction
#     #myreduce.recipename = 'testRegressionDark'
#     myreduce.runr()


import os
from pathlib import Path

from gempy.adlibrary import dataselect
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa : important to load adclass tags

from gempy.utils import logutils
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

    arms = ['BLUE', 'RED']
    exptimes = ["60s", "120s", "300s", "600s", "900s", "1200s", "1800s"]

    for arm in arms:
        for exptime in exptimes:
            only_darks = dataselect.select_data(all_files, tags=['RAW', 'DARK', arm, exptime])
        
            # Run reduce on all selected files
            myreduce = Reduce()
            myreduce.files.extend(only_darks)
            myreduce.drpkg = 'maroonxdr'
            # coment out this line for default reduction
            #myreduce.recipename = 'testRegressionDark'
            myreduce.runr()

    # Restore original working directory
    os.chdir(original_dir)

def test_dark_coefficients():

    # Get all files in the science_dir.  Change the path here to suit your installation.
    test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
    science_dir = test_path / 'science_dir'
    masterdark_path = science_dir / 'calibrations' / 'processed_dark'

    # Get all files in the science_dir
    all_files = list(Path(masterdark_path).glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    # Change working directory to science_dir
    original_dir = os.getcwd()
    os.chdir(science_dir)

    arms = ['BLUE', 'RED']

    for arm in arms:
        
        only_darks = dataselect.select_data(all_files, tags=['PROCESSED', 'DARK', arm], xtags=['DARK_COEFF']) 
    
        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_darks)
        myreduce.drpkg = 'maroonxdr'
        # coment out this line for default reduction
        myreduce.recipename = 'makeDarkCoefficients'
        myreduce.runr()

    # Restore original working directory
    os.chdir(original_dir)    


if __name__ == '__main__':

    test_reduce_dark()

    test_dark_coefficients()