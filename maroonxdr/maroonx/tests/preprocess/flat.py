"""Build master flats and measure blaze functions.

Reads debundled flat frames from ``$DRAGONS_TEST/preprocessed_files/`` (produced
by preprocess/bundle.py) and runs the ``makeProcessedFlatDFFFF`` recipe per arm:
stream separation (DFFFD vs DDDDF), stacking, straylight removal, stripe
finding/identification/definition, producing the DFFFF master flat consumed by
the wavecal and science staging steps. Then measures the blaze function on each
master flat.

Outputs are written under ``$DRAGONS_TEST/preprocessed_files/`` (and its
``calibrations/`` store).

Usage:
    python -m maroonxdr.maroonx.tests.preprocess.flat
"""

import os
from pathlib import Path

from gempy.adlibrary import dataselect
from gempy.utils import logutils
from recipe_system.reduction.coreReduce import Reduce

import maroonx_instruments  # noqa - import is necessary for dataselect
from maroonxdr.maroonx.tests.test_utils import change_cwd_context

ARMS = ['BLUE', 'RED']


def _get_dragons_test():
    p = os.environ.get('DRAGONS_TEST')
    if p is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')
    return Path(p)


def make_masterflats():
    """Build the DFFFF master flat for each arm from the debundled flats."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='flat.log', stomp=False)
        log = logutils.get_logger('flat.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            only_flats = dataselect.select_data(
                all_files, tags=['RAW', 'FLAT', arm]
            )
            myreduce = Reduce()
            myreduce.files.extend(only_flats)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeProcessedFlatDFFFF'
            myreduce.runr()


def make_blaze():
    """Measure the blaze function on each arm's master flat."""
    dragons_test = _get_dragons_test()
    preprocessed_dir = dragons_test / 'preprocessed_files'

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(preprocessed_dir):
        logutils.config(file_name='flat.log', stomp=False)
        log = logutils.get_logger('flat.log')
        log.setLevel('DEBUG')

        for arm in ARMS:
            selected_mflats = dataselect.select_data(
                all_files, tags=['PROCESSED', 'FLAT', arm]
            )
            myreduce = Reduce()
            myreduce.files.extend(selected_mflats)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeBlaze'
            myreduce.runr()


if __name__ == '__main__':

    # stack debundled flats into DFFFF master flats
    make_masterflats()

    # measure the blaze function on the master flats
    make_blaze()
