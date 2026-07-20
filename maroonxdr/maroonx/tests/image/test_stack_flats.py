"""Unit tests for the stackFlats primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

FDDDF_FIBER_SETUP = ['Flat lamp', 'Dark', 'Dark', 'Dark', 'Flat lamp']
LEVELS = [100.0, 200.0, 300.0]


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
