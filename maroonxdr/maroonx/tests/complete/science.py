"""Test the science reduction for MAROON-X data.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
complete/bundle.py) and writes reduced science outputs back to the same directory.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a reduced science
file. End users should use this test to test their installation. Only if there is an
error in this test (or other "complete" tests) should they use the echelle and image
unit tests to test their installation.

Make sure that you have created the darks and flats first (see dark.py and flat.py).
"""

import os
import shutil
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for astrodata
from maroonxdr.maroonx.tests.test_utils import change_cwd_context


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def complete_synthetic_darks_reduction():
    """Test reduction of synthetic darks for science frames."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    # Read debundled files from preprocessed_files/
    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_synth_dark.log', mode='debug', stomp=True)

        for arm in ['BLUE', 'RED']:
            selected_sci = dataselect.select_data(all_files, tags=['RAW', 'SCI', arm])

            # Run reduce on all selected files
            myreduce = Reduce()
            myreduce.files.extend(selected_sci)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeSyntheticDark'
            myreduce.runr()


def complete_science_reduction():
    """Test reduction of science frames across both arms (300s)."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    # Read debundled files from preprocessed_files/
    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_reduction.log', mode='debug', stomp=True)
        # logutils.config(file_name="test_science.log", stomp=True)
        # log = logutils.get_logger("test_science.log")
        # log.setLevel("DEBUG")

        for arm in ['BLUE', 'RED']:
            only_science = dataselect.select_data(
                all_files, tags=['RAW', 'SCI', arm, '300s']
            )

            myreduce = Reduce()
            myreduce.files.extend(only_science)
            myreduce.drpkg = 'maroonxdr'
            myreduce.runr()


def populate_inputs():
    """Copy science outputs from preprocessed_files/ to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # TODO: populate synth dark inputs

    # echelle_extraction/test_extraction
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_extraction' / 'inputs',
        [
            '20241124T041907Z_SOOOE_r_0300_reduced.fits',
            '20241124T041907Z_SOOOE_b_0300_reduced.fits',
        ],
    )

    # echelle_extraction/test_stripe_retrieval
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_stripe_retrieval' / 'inputs',
        [
            '20241124T041907Z_SOOOE_r_0300_test_stripes.fits',
            '20241124T041907Z_SOOOE_b_0300_test_stripes.fits',
        ],
    )


def _copy_files(src_dir, dst_dir, filenames):
    """Copy specific files from src_dir to dst_dir, creating dst_dir if needed."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    for f in filenames:
        src_file = src_dir / f
        if src_file.exists():
            shutil.copy2(src_file, dst_dir / f)
            print(f'  Copied {f} -> {dst_dir}')
        else:
            print(f'  WARNING: {src_file} not found, skipping')


if __name__ == '__main__':
    import sys

    complete_synthetic_darks_reduction()
    complete_science_reduction()

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs()
