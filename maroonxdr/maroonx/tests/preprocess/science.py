"""Build synthetic darks and reduce the science frames.

Reads debundled HD3651 science (SOOOE) frames from
``$DRAGONS_TEST/preprocessed_files/`` (produced by preprocess/bundle.py) and
runs two steps per arm:

1. ``makeSyntheticDark``: builds a synthetic dark matching the science exposure
   time from the dark coefficients (passed explicitly to ``createSyntheticDark``).
2. The default SCI recipe (``reduce``): stripe extraction, optimal extraction,
   etalon peak fitting, wavelength solution and drift correction, fiber
   combination, and barycentric correction.

All calibrations are passed explicitly per arm (hardcoded filenames read from
the working-directory copies in ``preprocessed_files/``), so the run does not
depend on caldb association: the master flat and synthetic dark go to
``extractStripes``, the dynamic wavecal to ``applyWavelengthSolution``, and the
dark coefficients to ``createSyntheticDark``.

Make sure that you have created the darks, flats, and wavecals first (see
dark.py, flat.py, and wavecal.py).

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.science
"""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for dataselect
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

ARMS = ['BLUE', 'RED']

# Calibrations produced by the upstream staging steps, per arm.
MASTERFLATS = {
    'BLUE': '20250701T171553Z_DDDDF_b_0007_DFFFF_flat.fits',
    'RED': '20250701T171553Z_DDDDF_r_0002_DFFFF_flat.fits',
}
DARK_COEFFS = {
    'BLUE': '20250707T164838Z_DDDDE_b_0120_darkCoefficients.fits',
    'RED': '20250707T164838Z_DDDDE_r_0120_darkCoefficients.fits',
}
SYNTH_DARKS = {
    'BLUE': '20250717T144308Z_SOOOE_b_0300_synth_dark.fits',
    'RED': '20250717T144308Z_SOOOE_r_0300_synth_dark.fits',
}
WAVECALS = {
    'BLUE': '20250717T163124Z_DEEEE_b_0010_wavecal.fits',
    'RED': '20250717T163124Z_DEEEE_r_0004_wavecal.fits',
}


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def make_synthetic_darks():
    """Build a synthetic dark for each arm's science exposure time."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='science.log', stomp=False)
        log = logutils.get_logger('science.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            selected_sci = dataselect.select_data(
                all_files, tags=['RAW', 'SCI', arm]
            )

            myreduce = Reduce()
            myreduce.files.extend(selected_sci)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeSyntheticDark'
            myreduce.uparms = {'createSyntheticDark:dark_coeff': DARK_COEFFS[arm]}
            myreduce.runr()


def make_science_reduction():
    """Reduce the HD3651 science frames for both arms."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='science.log', stomp=False)
        log = logutils.get_logger('science.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            only_science = dataselect.select_data(
                all_files, tags=['RAW', 'SCI', arm, '300s']
            )

            myreduce = Reduce()
            myreduce.files.extend(only_science)
            myreduce.drpkg = 'maroonxdr'
            myreduce.uparms = {
                'extractStripes:flat': MASTERFLATS[arm],
                'extractStripes:dark': SYNTH_DARKS[arm],
                'extractStripes:straylight_removal_fibers': [5],
                'applyWavelengthSolution:wavecal': WAVECALS[arm],
                'combineFibers:max_clips': 20,
            }
            myreduce.runr()


if __name__ == '__main__':

    # build synthetic darks for the science exposure times
    make_synthetic_darks()

    # reduce the science frames with all calibrations passed explicitly
    make_science_reduction()
