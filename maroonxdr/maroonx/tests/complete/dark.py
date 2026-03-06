"""Test the creation of master darks for MAROON-X data.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
complete/bundle.py) and writes master darks back to the same directory.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a calibrated dark
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.
"""

import itertools as it
import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.tests.test_utils import change_cwd_context


def _get_dragons_test():
    p = os.environ.get("DRAGONS_TEST")
    if p is None:
        raise RuntimeError("DRAGONS_TEST environment variable not set")
    return Path(p)


def complete_masterdark_reduction():
    """Test reduction of dark frames across all arms and exposure times."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    # Read debundled files from preprocessed_files/
    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name="test_dark.log", stomp=False)
        log = logutils.get_logger("test_dark.log")
        log.setLevel("DEBUG")

        arms = ['BLUE', 'RED']
        exptimes = ["60s", "120s", "300s", "600s", "900s", "1200s", "1800s"]

        for exptime, arm in it.product(exptimes, arms):

            only_darks = dataselect.select_data(
                all_files, tags=['RAW', 'DARK', arm, exptime])

            myreduce = Reduce()
            myreduce.files.extend(only_darks)
            myreduce.drpkg = 'maroonxdr'
            myreduce.runr()


def complete_dark_coeff_reduction():
    """Test creation of dark scaling coefficients from processed darks."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    # Re-scan to pick up newly produced master darks
    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name="test_dark.log", stomp=False)
        log = logutils.get_logger("test_dark.log")
        log.setLevel("DEBUG")

        for arm in ['BLUE', 'RED']:

            only_darks = dataselect.select_data(all_files,
                tags=['PROCESSED', 'DARK', arm], xtags=['DARK_COEFF', 'DARK_SYNTH'])

            myreduce = Reduce()
            myreduce.files.extend(only_darks)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeDarkCoefficients'
            myreduce.runr()


if __name__ == '__main__':

    complete_masterdark_reduction()

    complete_dark_coeff_reduction()
