"""
Script used to test the creation of 1D reduced spectra for MAROON-X wavecal files.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a 1D reduced spectra
file that can be used for dynamic wavelength calibration. End users should use this
test to test their installation.

To run this test, simply run the command "python wavecal_test.py" in the terminal, after
ensuring the correct path to science_dir in Path(). This test expects you to have a
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
def complete_wavecal_reduction():
    """Test reduction of wavelength calibration frames for both red and blue arms."""

    # Configure test logging
    logutils.config(file_name="test_wavecal.log", stomp=False)
    log = logutils.get_logger("test_wavecal.log")
    log.setLevel("DEBUG")

    # Get all files in the science_dir
    all_files = list(Path().glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()

    for arm in ['BLUE', 'RED']:
        only_wavecal = dataselect.select_data(all_files, tags=['RAW', 'WAVECAL', arm])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_wavecal)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()


if __name__ == '__main__':

    complete_wavecal_reduction()
