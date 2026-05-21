"""Run wavecal reduction on v2 (202507xx) etalon data.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
preprocess/bundle.py) and writes wavecal outputs back to the same directory.

Make sure that you have created the darks and flats first (see dark.py and flat.py).

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.wavecal [--populate-inputs] [--legacy-patch]
"""

import os
import shutil
import sys
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


def complete_wavecal_reduction(legacy_patch=False):
    """Reduce v2 etalon frames for both red and blue arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

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
            myreduce.uparms = {'extractStripes:legacy': legacy_patch}
            myreduce.runr()


def populate_inputs(legacy_patch=False):
    """Copy wavecal outputs to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # echelle_extraction/test_wavecal
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_wavecal' / 'inputs',
        [
            '20250717T163124Z_DEEEE_b_0010_wavecal.fits',
        ],
    )

    if not legacy_patch:
        return

    _copy_files(
        src,
        base / 'legacy_regression' / 'test_reduced_wavecal' / 'inputs',
        [
            '20250717T163124Z_DEEEE_b_0010_wavecal.fits',
            '20250717T163124Z_DEEEE_r_0004_wavecal.fits',
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

    legacy_patch = '--legacy-patch' in sys.argv[1:]
    complete_wavecal_reduction(legacy_patch=legacy_patch)

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs(legacy_patch=legacy_patch)
