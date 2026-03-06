"""Test bundle file reduction workflow.

Reads raw bundles from $DRAGONS_TEST/raw_files/ and writes debundled
single-arm files to $DRAGONS_TEST/preprocessed_files/.
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


def complete_bundle_reduction():
    """Test reduction of bundle FITS files containing both red and blue arms."""
    dragons_test = _get_dragons_test()
    raw_dir = dragons_test / 'raw_files'
    output_dir = dragons_test / 'preprocessed_files'
    output_dir.mkdir(exist_ok=True)

    # Read bundles from raw_files/
    all_files = sorted(str(p) for p in raw_dir.glob('*.fits'))
    only_bundles = dataselect.select_data(all_files, tags=['RAW', 'BUNDLE'])

    # Write debundled output to preprocessed_files/
    with change_cwd_context(output_dir):
        logutils.config(file_name='test_bundle.log', stomp=False)
        log = logutils.get_logger('test_bundle.log')
        log.setLevel('DEBUG')

        myreduce = Reduce()
        myreduce.files.extend(only_bundles)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()


def populate_inputs():
    """Copy debundled files from preprocessed_files/ to test inputs/ directories."""
    dragons_test = _get_dragons_test()
    src = dragons_test / 'preprocessed_files'
    base = dragons_test / 'maroonxdr' / 'maroonx'

    # bundle tests need raw bundles from raw_files/, not debundled — skip here
    # (bundle test create_inputs() downloads from archive directly)

    # image/test_file_sorting
    _copy_files(
        src,
        base / 'image' / 'test_file_sorting' / 'inputs',
        [
            '20241114T181028Z_DFFFD_r_0002.fits',
            '20241114T181028Z_DFFFD_b_0008.fits',
            '20241114T181815Z_DFFFD_b_0008.fits',
            '20241114T191006Z_DDDDF_b_0007.fits',
        ],
    )

    # image/test_image_orientation_corrector
    _copy_files(
        src,
        base / 'image' / 'test_image_orientation_corrector' / 'inputs',
        [
            '20241114T181815Z_DFFFD_r_0002.fits',
            '20241114T181959Z_DFFFD_b_0008.fits',
        ],
    )

    # image/test_ND_filter_check
    _copy_files(
        src,
        base / 'image' / 'test_ND_filter_check' / 'inputs',
        [
            '20241114T181028Z_DFFFD_r_0002.fits',
        ],
    )

    # image/test_var
    _copy_files(
        src,
        base / 'image' / 'test_var' / 'inputs',
        [
            '20241115T194624Z_DDDDE_r_0300.fits',
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

    complete_bundle_reduction()

    if '--populate-inputs' in sys.argv[1:]:
        populate_inputs()
