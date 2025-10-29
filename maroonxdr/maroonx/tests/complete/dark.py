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

import itertools as it
import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

# Get all files in the science_dir.
test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
science_dir = test_path / 'science_dir'


@change_cwd_context(science_dir)
def complete_masterdark_reduction():
    """Test reduction of dark frames across all arms and exposure times."""

    # Configure test logging
    logutils.config(file_name="test_dark.log", stomp=False)
    log = logutils.get_logger("test_dark.log")
    log.setLevel("DEBUG")

    # Get all files
    all_files = list(Path().glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    arms = ['BLUE', 'RED']
    exptimes = ["60s", "120s", "300s", "600s", "900s", "1200s", "1800s"]

    for exptime, arm in it.product(exptimes, arms):

        only_darks = dataselect.select_data(
            all_files, tags=['RAW', 'DARK', arm, exptime])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_darks)
        myreduce.drpkg = 'maroonxdr'
        # coment out this line for default reduction
        #myreduce.recipename = 'testRegressionDark'
        myreduce.runr()


@change_cwd_context(science_dir)
def complete_dark_coeff_reduction():
    """Test creation of dark scaling coefficients from processed darks."""

    # Configure test logging
    logutils.config(file_name="test_dark.log", stomp=False)
    log = logutils.get_logger("test_dark.log")
    log.setLevel("DEBUG")

    # Get all files
    masterdark_path = Path() / 'calibrations' / 'processed_dark'
    all_files = list(Path(masterdark_path).glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    arms = ['BLUE', 'RED']

    for arm in arms:

        only_darks = dataselect.select_data(all_files, 
            tags=['PROCESSED', 'DARK', arm], xtags=['DARK_COEFF', 'DARK_SYNTH'])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_darks)
        myreduce.drpkg = 'maroonxdr'
        # coment out this line for default reduction
        myreduce.recipename = 'makeDarkCoefficients'
        myreduce.runr()


if __name__ == '__main__':

    complete_masterdark_reduction()

    complete_dark_coeff_reduction()
