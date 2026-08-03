"""Regression tests for the removeStrayLight primitive.

The input frames are written into inputs/ by ``create_inputs_recipe()``, which
runs the makeStrayLightCheck recipe on the debundled flats: stacked flat
streams staged just before straylight removal, with the blessed result
stored as the STRAYLIGHT_DIFFERENCE extension.

Usage:
    python -m maroonxdr.maroonx.tests.image.test_remove_straylight --create-inputs
"""

import os

import numpy as np
import pytest
from gempy.utils import logutils

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

# Snapshot frames written by the makeStrayLightCheck recipe, per arm and stream.
datasets = [
    '20250701T171553Z_DDDDF_b_0007_straylight_flat.fits',
    '20250701T171553Z_DDDDF_r_0002_straylight_flat.fits',
    '20250701T170101Z_DFFFD_b_0008_straylight_flat.fits',
    '20250701T170101Z_DFFFD_r_0002_straylight_flat.fits',
]


# -- Tests ---------------------------------------------------------------------
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filename', datasets)
def test_remove_straylight(filename, path_to_inputs, change_working_dir):
    """removeStrayLight output is stable against the blessed snapshot."""
    ad = astrodata.open(os.path.join(path_to_inputs, filename))
    expected = ad[0].data + ad[0].STRAYLIGHT_DIFFERENCE

    with change_working_dir():
        logutils.config(file_name=f'log_{filename.replace(".fits", "")}.txt')
        ad_out = MAROONX([]).removeStrayLight([ad], report=False).pop()

    # dtype is float32
    np.testing.assert_allclose(ad_out[0].data, expected, rtol=1e-6, atol=1e-4)

    # sanity check: the inter-order background is removed
    inter_order = ad_out[0].INDEX_FIBER == 0
    assert np.abs(np.median(ad_out[0].data[inter_order])) < 1.0


# -- Recipe to create the input snapshots --------------------------------------
def create_inputs_recipe():
    """Stage the straylight snapshot frames into inputs/."""
    from pathlib import Path
    from gempy.adlibrary import dataselect
    from recipe_system.reduction.coreReduce import Reduce
    from maroonxdr.maroonx.tests.test_utils import change_cwd_context

    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST not set')

    preprocessed_dir = Path(root) / 'preprocessed_files'
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    inputs_dir = (
        Path(root) / 'maroonxdr' / 'maroonx' / 'image' / module_name / 'inputs'
    )
    inputs_dir.mkdir(parents=True, exist_ok=True)

    all_files = sorted(str(p) for p in preprocessed_dir.glob('*.fits'))

    with change_cwd_context(inputs_dir):
        logutils.config(file_name='create_inputs.log', stomp=False)

        for arm in ['BLUE', 'RED']:
            only_flats = dataselect.select_data(
                all_files, tags=['RAW', 'FLAT', arm]
            )
            myreduce = Reduce()
            myreduce.files.extend(only_flats)
            myreduce.drpkg = 'maroonxdr'
            myreduce.recipename = 'makeStrayLightCheck'
            myreduce.runr()


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
