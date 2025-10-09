"""Test bundle file reduction workflow."""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

logutils.config(file_name="test_reduction.log", stomp=False)
log = logutils.get_logger("test_reduction.log")
log.setLevel("DEBUG")


def test_reduce_bundle():
    """Test reduction of bundle FITS files containing both red and blue arms."""
    # Get all files in the science_dir.  Change the path here to suit your installation.
    test_path = Path(os.environ.get("MAROONX_DRAGONS_TEST"))
    science_dir = test_path / 'science_dir'

    # Get all files in the science_dir
    all_files = list(Path(science_dir).glob('*.fits'))
    all_files = [str(p) for p in all_files]
    all_files.sort()
    only_bundles = dataselect.select_data(all_files, tags=['BUNDLE'])

    # Change working directory to science_dir
    original_dir = Path.cwd()
    os.chdir(science_dir)

    # Run reduce on all selected files
    myreduce = Reduce()
    myreduce.files.extend(only_bundles)
    myreduce.drpkg = 'maroonxdr'
    myreduce.runr()

    # Restore original working directory
    os.chdir(original_dir)


if __name__ == '__main__':

    test_reduce_bundle()
