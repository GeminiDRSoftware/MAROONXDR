"""
Script used to test the creation of flats for MAROON-X data.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a calibrated flat
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.

To run this test, simply run the command "python flat_test.py" in the terminal, after
ensuring the correct path to science_dir in Path().  This test expects you to have a
science_dir in the root directory of the installation (you have to make this directory)
and to have got the test fits files from Kathleen.
"""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.tests.test_utils import change_cwd_context


test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
science_dir = test_path / 'science_dir'


@change_cwd_context(science_dir)
def complete_masterflat_reduction():
    """Test reduction of flat frames for both red and blue arms."""

    # Configure test logging
    logutils.config(file_name="test_flat.log", stomp=False)
    log = logutils.get_logger("test_flat.log")
    log.setLevel("DEBUG")

    # Get all files in the science_dir
    all_files = list(Path().glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    for arm in ['BLUE', 'RED']:
        only_flats = dataselect.select_data(all_files, tags=['RAW', 'FLAT', arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_flats)
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeProcessedFlatDFFFF'
        myreduce.runr()


if __name__ == '__main__':

    complete_masterflat_reduction()
