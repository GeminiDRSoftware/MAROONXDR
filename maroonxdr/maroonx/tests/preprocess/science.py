"""Run synthetic dark creation and science reduction on v2 (202507xx) data.

Reads debundled files from $DRAGONS_TEST/preprocessed_files/ (produced by
preprocess/bundle.py), creates synthetic darks, and runs the full science
reduction.

Make sure that you have created the darks, flats, and wavecal first
(see dark.py, flat.py, and wavecal.py).

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.science [--populate-inputs] [--legacy-patch]
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


def complete_synthetic_darks_reduction():
    """Create synthetic darks for v2 science frames."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_synth_dark.log', mode='debug', stomp=True)

        for arm in ['BLUE', 'RED']:
            selected_sci = dataselect.select_data(all_files, tags=['RAW', 'SCI', arm])

            myreduce = Reduce()
            myreduce.files.extend(selected_sci)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeSyntheticDark'
            myreduce.runr()


def complete_science_reduction(legacy_patch=False):
    """Reduce v2 science frames for both red and blue arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_science.log', mode='debug', stomp=True)

        for arm in ['BLUE', 'RED']:
            only_science = dataselect.select_data(
                all_files, tags=['RAW', 'SCI', arm, '300s']
            )

            myreduce = Reduce()
            myreduce.files.extend(only_science)
            myreduce.drpkg = 'maroonxdr'
            myreduce.uparms = {
                'extractStripes:legacy': legacy_patch,
                'combineFibers:max_clips': 20,
                'extractStripes:straylight_removal_fibers': [5],
            }
            myreduce.runr()


def complete_stripe_extraction_check():
    """Create test_stripes files by running makeStripeExtractionCheck recipe."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='test_stripe_check.log', mode='debug', stomp=True)

        for arm in ['BLUE', 'RED']:
            only_science = dataselect.select_data(
                all_files, tags=['RAW', 'SCI', arm, '300s']
            )

            myreduce = Reduce()
            myreduce.files.extend(only_science)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeStripeExtractionCheck'
            myreduce.runr()


def populate_inputs(legacy_patch=False):
    """Copy synth dark and reduced science outputs to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # echelle_extraction/test_extraction: needs reduced science (both arms)
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_extraction' / 'inputs',
        [
            '20250717T144308Z_SOOOE_b_0300_reduced.fits',
            '20250717T144308Z_SOOOE_r_0300_reduced.fits',
        ],
    )

    # echelle_extraction/test_stripe_retrieval: needs test_stripes files
    _copy_files(
        src,
        base / 'echelle_extraction' / 'test_stripe_retrieval' / 'inputs',
        [
            '20250717T144308Z_SOOOE_b_0300_test_stripes.fits',
            '20250717T144308Z_SOOOE_r_0300_test_stripes.fits',
        ],
    )

    if not legacy_patch:
        return

    # legacy_regression/test_masterdark: needs synth darks
    dark_src = src / 'calibrations' / 'processed_dark'
    _copy_files(
        dark_src,
        base / 'legacy_regression' / 'test_masterdark' / 'inputs',
        [
            '20250717T144308Z_SOOOE_b_0300_synth_dark.fits',
            '20250717T144308Z_SOOOE_r_0300_synth_dark.fits',
        ],
    )

    # legacy_regression/test_reduced_science: needs reduced science (both arms)
    _copy_files(
        src,
        base / 'legacy_regression' / 'test_reduced_science' / 'inputs',
        [
            '20250717T144308Z_SOOOE_b_0300_reduced.fits',
            '20250717T144308Z_SOOOE_r_0300_reduced.fits',
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

    complete_synthetic_darks_reduction()
    complete_science_reduction(legacy_patch=legacy_patch)
    complete_stripe_extraction_check()

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs(legacy_patch=legacy_patch)
