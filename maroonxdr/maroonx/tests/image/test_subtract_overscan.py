"""Unit tests for the subtractOverscan primitive."""

import numpy as np
import pytest

from maroonxdr.maroonx.primitives_maroonx_2D import MAROONX

from . import make_frame

# The overscan sections are read from a static lookup, so the frame needs the
# real detector size for them to address the right pixels.
DETECTOR_SHAPE = (4400, 4400)
SIGNAL = 100.0
BIAS = 11.0


# -- Fixtures ------------------------------------------------------------------
@pytest.fixture(params=['RED', 'BLUE'])
def ad_biased(request):
    """Illuminated pixels carry SIGNAL on top of BIAS; the overscan carries BIAS only."""
    array = np.full(DETECTOR_SHAPE, SIGNAL + BIAS, dtype=np.float32)
    ad = make_frame(request.param, array)
    for sec in ad.subtract_overscan_section()[0]:
        ad[0].data[sec.asslice()] = BIAS
    return ad


# -- Tests ---------------------------------------------------------------------
def test_subtractOverscan(ad_biased):
    """Overscan mean is subtracted from each array quadrant."""
    ad = ad_biased
    osec = ad.subtract_overscan_section()[0]

    result = MAROONX([]).subtractOverscan([ad])[0]
    data = result[0].data

    for osec_i in osec:
        np.testing.assert_allclose(data[osec_i.asslice()], 0.0, atol=1e-5)


def test_subtractOverscan_removes_bias_from_data(ad_biased):
    """The overscan bias is also removed from the illuminated pixels.

    The overscan carries BIAS, so that is what is subtracted from each quadrant,
    leaving the illuminated pixels at SIGNAL.
    """
    ad = ad_biased
    osec = ad.subtract_overscan_section()[0]

    # Everything the overscan strips do not cover
    is_data = np.ones(DETECTOR_SHAPE, dtype=bool)
    for osec_i in osec:
        is_data[osec_i.asslice()] = False

    result = MAROONX([]).subtractOverscan([ad])[0]
    data = result[0].data

    np.testing.assert_allclose(data[is_data], SIGNAL, rtol=1e-6)
