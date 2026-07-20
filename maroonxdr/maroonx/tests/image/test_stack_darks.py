"""Unit tests for the stackDarks primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

LEVELS = [100.0, 200.0, 300.0]


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture(params=['RED', 'BLUE'])
def ad_darks(request):
    """Three DDDDE darks with a different uniform flux level each."""
    ad_list = []
    for level in LEVELS:
        ad_list.append(
            make_frame(request.param, np.full((32, 32), level, dtype=np.float32))
        )
    return ad_list

# -- Tests ---------------------------------------------------------------------
def test_stackDarks(ad_darks):
    """Frames with different flux levels are scaled to the first frame before combining."""
    result = MAROONX([]).stackDarks(ad_darks, scale_mode='first_frame')[0]

    np.testing.assert_allclose(result[0].data, LEVELS[0], rtol=1e-6)
    assert result.phu.get('NCOMBINE') == len(LEVELS)


def test_stackDarks_mixed_exposure_times(ad_darks):
    """Darks of unequal exposure time cannot be combined."""
    # the first frame sets the reference exptime at 300
    # test when the last file is wrong
    ad_darks[-1].phu['EXPTIME'] = 600.0

    with pytest.raises(ValueError, match='not of equal exposure time'):
        MAROONX([]).stackDarks(ad_darks, scale_mode='first_frame')
