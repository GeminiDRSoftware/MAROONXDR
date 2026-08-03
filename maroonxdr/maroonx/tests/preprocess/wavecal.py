"""Build dynamic wavelength calibrations from etalon frames.

Reads debundled etalon (DEEEE) frames from ``$DRAGONS_TEST/preprocessed_files/``
(produced by preprocess/bundle.py) and runs the default WAVECAL recipe
(``makeDynamicWavecal``) per arm: stripe extraction, box extraction, etalon
peak fitting, static solution lookup, and the dynamic etalon fit.

The master flat produced by preprocess/flat.py is passed explicitly to
``extractStripes`` (hardcoded per arm in ``MASTERFLATS``, read from the
working-directory copy in ``preprocessed_files/``), so the run does not depend
on caldb association for the flat. Dark subtraction is skipped for all fibers
by default in this recipe, so no dark is used.

Make sure that you have created the darks and flats first (see dark.py and
flat.py).

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.wavecal
"""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for dataselect
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

ARMS = ['BLUE', 'RED']

# Master flats produced by preprocess/flat.py, per arm.
MASTERFLATS = {
    'BLUE': '20250701T171553Z_DDDDF_b_0007_DFFFF_flat.fits',
    'RED': '20250701T171553Z_DDDDF_r_0002_DFFFF_flat.fits',
}


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def make_wavecals():
    """Build the dynamic wavelength calibration for each arm's etalon frame."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='wavecal.log', stomp=False)
        log = logutils.get_logger('wavecal.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            only_wavecal = dataselect.select_data(
                all_files, tags=['RAW', 'WAVECAL', arm]
            )
            myreduce = Reduce()
            myreduce.files.extend(only_wavecal)
            myreduce.drpkg = 'maroonxdr'
            myreduce.uparms = {'extractStripes:flat': MASTERFLATS[arm]}
            myreduce.runr()


if __name__ == '__main__':

    # reduce etalons per arm
    make_wavecals()
