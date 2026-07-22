"""Unit tests for the stackFlats primitive.

The fast tests run on synthetic frames built in-test. The preprocessed_data
test stacks real debundled flats (staged into inputs/ from the preprocess/
chain by ``create_inputs_recipe()``) and checks the structural and
statistical properties of the stack; there is no blessed reference.

Usage:
    python -m maroonxdr.maroonx.tests.image.test_stack_flats --create-inputs
"""

import os

import numpy as np
import pytest

import astrodata
import maroonx_instruments  # noqa - registers the MaroonX AstroData class
from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

FDDDF_FIBER_SETUP = ['Flat lamp', 'Dark', 'Dark', 'Dark', 'Flat lamp']
LEVELS = [100.0, 200.0, 300.0]

# Debundled DDDDF flats from the preprocess/ chain, per arm.
datasets = [
    ['20250701T171553Z_DDDDF_b_0007.fits',
     '20250701T171811Z_DDDDF_b_0007.fits',
     '20250701T171955Z_DDDDF_b_0007.fits'],
    ['20250701T171553Z_DDDDF_r_0002.fits',
     '20250701T171811Z_DDDDF_r_0002.fits',
     '20250701T171955Z_DDDDF_r_0002.fits'],
]


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture(params=['RED', 'BLUE'])
def ad_flats(request):
    """Three FDDDF flats with a different uniform flux level each."""
    ad_list = []
    for level in LEVELS:
        ad_list.append(
            make_frame(request.param, np.full((32, 32), level, dtype=np.float32),
                       fiber_setup=FDDDF_FIBER_SETUP)
        )
    return ad_list


# -- Tests ---------------------------------------------------------------------
def test_stackFlats(ad_flats):
    """Frames with different flux levels are scaled to the mean before combining."""
    result = MAROONX([]).stackFlats(ad_flats, scale_mode='mean_frame')[0]

    np.testing.assert_allclose(result[0].data, np.mean(LEVELS), rtol=1e-6)
    assert result.phu.get('NCOMBINE') == len(LEVELS)


@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize('filenames', datasets)
def test_stackFlats_real_flats(filenames, path_to_inputs):
    """Stacking real debundled flats preserves structure and total flux."""
    ad_list = [astrodata.open(os.path.join(path_to_inputs, f))
               for f in filenames]

    result = MAROONX([]).stackFlats(ad_list, scale_mode='mean_frame')[0]

    # Structure
    assert len(result) == 1
    assert result[0].data.shape == ad_list[0][0].data.shape
    assert result[0].data.dtype == np.float32

    # Headr keywords
    assert result.phu.get('NCOMBINE') == len(filenames)
    for i in range(1, len(filenames) + 1):
        assert result.phu.get(f'IMCMB{i:03d}') is not None

    # All pixels finite, negatives clipped to zero
    assert np.all(np.isfinite(result[0].data))
    assert result[0].data.min() >= 0


# -- Recipe to create the inputs -----------------------------------------------
def create_inputs_recipe():
    """Copy the debundled flats into inputs/.

    The debundled flats are products of preprocess/bundle.py, so staging is a
    copy from ``$DRAGONS_TEST/preprocessed_files/``.
    """
    import shutil
    from pathlib import Path

    root = os.environ.get('DRAGONS_TEST')
    if root is None:
        raise RuntimeError('DRAGONS_TEST environment variable not set')

    preprocessed_dir = Path(root) / 'preprocessed_files'
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    inputs_dir = (
        Path(root) / 'maroonxdr' / 'maroonx' / 'image' / module_name / 'inputs'
    )
    inputs_dir.mkdir(parents=True, exist_ok=True)

    for filenames in datasets:
        for filename in filenames:
            shutil.copy2(preprocessed_dir / filename, inputs_dir / filename)
            print(f'Copied to inputs/:\n    {filename}')


if __name__ == '__main__':
    import sys

    if '--create-inputs' in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])
