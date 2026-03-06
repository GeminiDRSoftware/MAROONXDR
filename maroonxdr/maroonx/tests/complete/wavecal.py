"""Test the creation of 1D reduced spectra for MAROON-X wavecal files.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
complete/bundle.py) and writes wavecal outputs back to the same directory.

It does not rely on pytest, and does not produce a success or fail output like pytest
does. Instead, if the reduce runs successfully, this will produce a 1D reduced spectra
file that can be used for dynamic wavelength calibration. End users should use this
test to test their installation.
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


def complete_wavecal_reduction():
    """Test reduction of wavelength calibration frames for both red and blue arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    # Read debundled files from preprocessed_files/
    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_wavecal.log', stomp=False)
        log = logutils.get_logger('test_wavecal.log')
        log.setLevel('DEBUG')

        for arm in ['BLUE', 'RED']:
            only_wavecal = dataselect.select_data(
                all_files, tags=['RAW', 'WAVECAL', arm]
            )

            myreduce = Reduce()
            myreduce.files.extend(only_wavecal)
            myreduce.drpkg = 'maroonxdr'
            myreduce.runr()


def populate_inputs():
    """Copy wavecal outputs from preprocessed_files/ to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # echelle_extraction/test_wavecal
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_wavecal' / 'inputs',
        [
            '20241124T030227Z_DEEEE_b_0030_wavecal.fits',
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

    complete_wavecal_reduction()

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs()
