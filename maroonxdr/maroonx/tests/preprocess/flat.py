"""Run master flat reduction on v2 (202507xx) data.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
preprocess/bundle.py) and writes master flats back to the same directory.

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.flat [--populate-inputs] [--legacy-patch]
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


def complete_masterflat_reduction(legacy_patch=False):
    """Test reduction of flat frames for both red and blue arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    # Legacy glob() ordering from masterflat FITS HISTORY headers.
    # DDDDF files first, then DFFFD — separateFlatStreams preserves
    # within-group order.
    _legacy_order = {
        'BLUE': [
            '20250701T172509Z_DDDDF_b_0007.fits',
            '20250701T172324Z_DDDDF_b_0007.fits',
            '20250701T172140Z_DDDDF_b_0007.fits',
            '20250701T171955Z_DDDDF_b_0007.fits',
            '20250701T171811Z_DDDDF_b_0007.fits',
            '20250701T171553Z_DDDDF_b_0007.fits',
            '20250701T170537Z_DFFFD_b_0008.fits',
            '20250701T170101Z_DFFFD_b_0008.fits',
            '20250701T171051Z_DFFFD_b_0008.fits',
            '20250701T170906Z_DFFFD_b_0008.fits',
            '20250701T170353Z_DFFFD_b_0008.fits',
            '20250701T170721Z_DFFFD_b_0008.fits',
        ],
        'RED': [
            '20250701T171955Z_DDDDF_r_0002.fits',
            '20250701T172140Z_DDDDF_r_0002.fits',
            '20250701T171553Z_DDDDF_r_0002.fits',
            '20250701T172324Z_DDDDF_r_0002.fits',
            '20250701T171811Z_DDDDF_r_0002.fits',
            '20250701T172509Z_DDDDF_r_0002.fits',
            '20250701T170906Z_DFFFD_r_0002.fits',
            '20250701T170721Z_DFFFD_r_0002.fits',
            '20250701T170101Z_DFFFD_r_0002.fits',
            '20250701T170353Z_DFFFD_r_0002.fits',
            '20250701T171051Z_DFFFD_r_0002.fits',
            '20250701T170537Z_DFFFD_r_0002.fits',
        ],
    }

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_flat.log', stomp=False)
        log = logutils.get_logger('test_flat.log')
        log.setLevel('DEBUG')

        for arm in ['BLUE', 'RED']:
            if legacy_patch:
                only_flats = [
                    str(preprocessed_dir / f) for f in _legacy_order[arm]
                ]
            else:
                only_flats = dataselect.select_data(
                    all_files, tags=['RAW', 'FLAT', arm]
                )

            myreduce = Reduce()
            myreduce.files.extend(only_flats)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeProcessedFlatDFFFF'
            myreduce.uparms = {'removeStrayLight:legacy': legacy_patch}
            myreduce.runr()


def complete_blaze_reduction():
    """Measure Blaze function on master flats for both red and blue arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_flat.log', stomp=False)
        log = logutils.get_logger('test_flat.log')
        log.setLevel('DEBUG')

        for arm in ['BLUE', 'RED']:
            selected_mflats = dataselect.select_data(
                all_files, tags=['PROCESSED', 'FLAT', arm]
            )
            myreduce = Reduce()
            myreduce.files.extend(selected_mflats)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'measureBlaze'
            myreduce.runr()


def complete_straylight_prep():
    """Prepare straylight-check flats for the straylight removal unit test."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    fdddf_file = str(preprocessed_dir / '20250701T171553Z_DDDDF_b_0007.fits')
    dfffd_file = str(preprocessed_dir / '20250701T170101Z_DFFFD_b_0008.fits')

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_flat.log', stomp=False)
        log = logutils.get_logger('test_flat.log')
        log.setLevel('DEBUG')

        myreduce = Reduce()
        myreduce.files.extend([dfffd_file, fdddf_file])
        myreduce.drpkg = 'maroonxdr'
        myreduce.recipename = 'makeStrayLightCheck'
        myreduce.runr()


def populate_inputs(legacy_patch=False):
    """Copy flat outputs from preprocessed_files/ to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    cal_src = src / 'calibrations' / 'processed_flat'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # image/test_stripe_finding: needs the red DFFFF master flat
    _copy_files(
        cal_src,
        base / 'image' / 'test_stripe_finding' / 'inputs',
        ['20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits'],
    )

    # image/test_stray_light_removal: needs straylight flats
    _copy_files(
        src,
        base / 'image' / 'test_stray_light_removal' / 'inputs',
        [
            '20250701T170101Z_DFFFD_b_0008_straylight_flat.fits',
            '20250701T171553Z_DDDDF_b_0007_straylight_flat.fits',
        ],
    )

    # echelle_extraction/test_measure_blaze: needs both arm DFFFF master flats
    _copy_files(
        cal_src,
        base / 'echelle_extraction' / 'test_measure_blaze' / 'inputs',
        [
            '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
        ],
    )

    if not legacy_patch:
        return

    # legacy_regression/test_masterflat: needs DFFFF processed flats (both arms)
    _copy_files(
        cal_src,
        base / 'legacy_regression' / 'test_masterflat' / 'inputs',
        [
            '20250701T172509Z_DDDDF_b_0007_DFFFF_flat.fits',
            '20250701T171955Z_DDDDF_r_0002_DFFFF_flat.fits',
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
    complete_masterflat_reduction(legacy_patch=legacy_patch)

    complete_blaze_reduction()

    complete_straylight_prep()

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs(legacy_patch=legacy_patch)
