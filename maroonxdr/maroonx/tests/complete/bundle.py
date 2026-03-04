"""Test bundle file reduction workflow."""

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


def complete_bundle_reduction():
    """Test reduction of bundle FITS files containing both red and blue arms."""
    science_dir = _get_science_dir()

    with change_cwd_context(science_dir):
        # Configure test logging
        logutils.config(file_name="test_bundle.log", stomp=False)
        log = logutils.get_logger("test_bundle.log")
        log.setLevel("DEBUG")

        # Get all files in the science_dir
        all_files = list(Path().glob('*.fits'))
        all_files = [str(p) for p in all_files]
        all_files.sort()
        only_bundles = dataselect.select_data(all_files, tags=['RAW', 'BUNDLE'])

        # Run reduce on all selected files
        myreduce = Reduce()
        myreduce.files.extend(only_bundles)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()


if __name__ == '__main__':

    complete_bundle_reduction()
