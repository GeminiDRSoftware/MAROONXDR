"""Fetch raw MaroonX bundles from the Gemini Archive and debundle them.

Downloads the raw two-arm bundle frames listed in ``MANIFEST`` from the Gemini
Archive into ``$DRAGONS_TEST/raw_files/`` (cached; skipped if already present),
then runs the default BUNDLE recipe (``splitBundle``) on each to produce
single-arm debundled frames in ``$DRAGONS_TEST/preprocessed_files/``.

The debundled frames are the shared inputs consumed by the dark, flat,
wavecal, and science staging steps.

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.bundle
"""

import os
import warnings
from pathlib import Path
from urllib.error import HTTPError

from astrodata.testing import download_from_archive
from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for dataselect
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

# Raw two-arm bundle frames on the Gemini Archive, grouped by calibration kind.
# Single source of truth for the raw file set.
MANIFEST = {
    'FLAT': [
        'N20250701M6126.fits',
        'N20250701M6143.fits',
        'N20250701M6154.fits',
        'N20250701M6164.fits',
        'N20250701M6175.fits',
        'N20250701M6185.fits',
        'N20250701M6215.fits',
        'N20250701M6229.fits',
        'N20250701M6240.fits',
        'N20250701M6250.fits',
        'N20250701M6260.fits',
        'N20250701M6271.fits',
    ],
    'DARK': [
        'N20250707M6052.fits',
        'N20250707M6074.fits',
        'N20250707M6096.fits',
        'N20250707M6119.fits',
        'N20250707M6141.fits',
        'N20250707M6164.fits',
        'N20250707M6180.fits',
        'N20250707M6197.fits',
        'N20250707M6213.fits',
        'N20250707M6230.fits',
        'N20250707M6246.fits',
        'N20250707M6287.fits',
        'N20250707M6327.fits',
        'N20250707M6368.fits',
        'N20250707M6408.fits',
        'N20250707M6449.fits',
        'N20250707M6520.fits',
        'N20250707M6590.fits',
        'N20250707M6660.fits',
        'N20250707M6731.fits',
        'N20250707M6802.fits',
        'N20250707M6902.fits',
        'N20250707M7002.fits',
        'N20250707M7102.fits',
        'N20250707M7203.fits',
        'N20250707M7304.fits',
        'N20250707M7434.fits',
        'N20250707M7564.fits',
        'N20250707M7695.fits',
        'N20250707M7826.fits',
        'N20250707M7956.fits',
        'N20250707M8147.fits',
        'N20250707M8337.fits',
        'N20250707M8527.fits',
        'N20250708M0078.fits',
    ],
    'WAVECAL': [
        'N20250717M5948.fits',
    ],
    'SCIENCE': [
        'N20250717M5299.fits',
    ],
}


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def download_raw_bundles():
    """Download every MANIFEST bundle into $DRAGONS_TEST/raw_files/ (cached).

    Returns
    -------
    dict
        Mapping of filename to local path, or None where the archive responded
        with HTTP 403 (proprietary data).
    """
    raw_dir = _get_dragons_test() / 'raw_files'

    paths = {}
    for filenames in MANIFEST.values():
        for filename in filenames:
            # Skip the archive query for cached files: the md5 check in
            # download_from_archive is broken for MaroonX filenames (GOA does
            # not filter its jsonfilelist response by them) and re-downloads
            # every file.
            local = raw_dir / filename
            if local.exists():
                paths[filename] = str(local)
                continue
            try:
                paths[filename] = download_from_archive(
                    filename, sub_path='raw_files', env_var='DRAGONS_TEST'
                )
            except HTTPError as e:
                warnings.warn(f'{filename}: {e}')
                paths[filename] = None
    return paths


def debundle():
    """Split raw bundles into single-arm frames in preprocessed_files/."""
    dragons_test = _get_dragons_test()
    raw_dir = dragons_test / 'raw_files'
    output_dir = dragons_test / 'preprocessed_files'
    output_dir.mkdir(exist_ok=True)

    all_files = sorted(str(p) for p in raw_dir.glob('N2025*.fits'))
    only_bundles = dataselect.select_data(all_files, tags=['RAW', 'BUNDLE'])

    with change_cwd_context(output_dir):
        logutils.config(file_name='bundle.log', stomp=False)
        log = logutils.get_logger('bundle.log')
        log.setLevel('DEBUG')

        myreduce = Reduce()
        myreduce.files.extend(only_bundles)
        myreduce.drpkg = 'maroonxdr'
        myreduce.runr()


if __name__ == '__main__':

    # download manifest into $DRAGONS_TEST/raw_files/
    download_raw_bundles()
    
    # debundle raw bundles into $DRAGONS_TEST/preprocessed_files/
    debundle()
