"""Build master darks and dark coefficients.

Reads debundled dark frames from ``$DRAGONS_TEST/preprocessed_files/`` (produced
by preprocess/bundle.py), stacks them into one master dark per (arm, exposure
time) via the default DARK recipe, then fits the per-pixel log-linear dark
coefficients from those master darks.

Outputs are written under ``$DRAGONS_TEST/preprocessed_files/`` (and its
``calibrations/`` store). The dark chain has no upstream calibration dependency,
so no calibrations are injected.

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.dark
"""

import itertools as it
import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for dataselect
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

# Dark exposure times present in the manifest (dataselect tags), and arms.
EXPTIMES = ['60s', '120s', '300s', '600s', '900s', '1200s', '1800s']
ARMS = ['BLUE', 'RED']


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def make_masterdarks():
    """Stack debundled darks into one master dark per (arm, exposure time)."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='dark.log', stomp=False)
        log = logutils.get_logger('dark.log')
        log.setLevel('DEBUG')

        for exptime, arm in it.product(EXPTIMES, ARMS):
            only_darks = dataselect.select_data(
                all_files, tags=['RAW', 'DARK', arm, exptime]
            )
            myreduce = Reduce()
            myreduce.files.extend(only_darks)
            myreduce.drpkg = 'maroonxdr'
            myreduce.runr()


def make_dark_coefficients():
    """Fit per-pixel log-linear dark coefficients from the master darks."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='dark.log', stomp=False)
        log = logutils.get_logger('dark.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            only_darks = dataselect.select_data(
                all_files,
                tags=['PROCESSED', 'DARK', arm],
                xtags=['DARK_COEFF', 'DARK_SYNTH'],
            )
            myreduce = Reduce()
            myreduce.files.extend(only_darks)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeDarkCoefficients'
            myreduce.runr()


if __name__ == '__main__':
    
    # create the master darks and store in caldb
    make_masterdarks()

    # create the dark coefficients and store in caldb
    make_dark_coefficients()
