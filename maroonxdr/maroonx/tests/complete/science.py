"""
Script used to test the science reduction for MAROON-X data.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a reduced science
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.

To run this test, simply run the command "python science_test.py" in the terminal.
This test expects you to have a science_dir in the root directory of the installation
(you have to make this directory) and to have got the test fits files from Kathleen.
You also have to make sure that you have created the darks and flats first. Make sure
these flats and darks are in a calibrations directory in the root directory of the
installation, and in processed_dark and processed_flat subdirectories, respectively.
If you need to change the paths, this can be done in primitives_maroonx_echelle.py.
"""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.tests.test_utils import change_cwd_context


def _get_science_dir():
    p = os.environ.get("DRAGONS_TEST")
    if p is None:
        raise RuntimeError("DRAGONS_TEST environment variable not set")
    return Path(p) / 'science_dir'


def complete_science_reduction():
    """Test reduction of science frames across both arms (300s)."""
    science_dir = _get_science_dir()

    with change_cwd_context(science_dir):
        # Configure test logging
        logutils.config(file_name="test_science.log", stomp=False)
        log = logutils.get_logger("test_science.log")
        log.setLevel("DEBUG")

        # Get all files
        all_files = list(Path().glob('*.fits'))
        all_files = [str(p) for p in all_files]
        all_files.sort()

        for arm in ['BLUE', 'RED']:
            only_science = dataselect.select_data(all_files,
                tags=['RAW', 'SCI', arm, '300s'])

            # Run reduce on all selected files
            myreduce = Reduce()
            myreduce.files.extend(only_science)
            myreduce.drpkg = 'maroonxdr'
            myreduce.runr()


if __name__ == '__main__':

    complete_science_reduction()
